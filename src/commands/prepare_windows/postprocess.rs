use crate::commands::prepare_windows::{
    config::{DedupKeep, DistanceTiesPolicy, MergeScope},
    prepare_windows::FinalWindow,
};
use fxhash::FxHashMap;
use std::cmp::Ordering;

/// Candidate retained while resolving spacing conflicts.
///
/// When enforcing a minimum distance the algorithm collects overlapping windows
/// and later chooses a winner according to the configured policy.
///
/// The candidate stores the index into the original window vector plus optional
/// score and length metadata used by tie-breaking strategies.
#[derive(Clone)]
struct Candidate {
    window_idx: usize,
    score: Option<f32>,
    length: u32,
}

/// Deduplicate identical windows according to the configured policy.
///
/// After sorting we can eliminate exact duplicates, keeping either the first,
/// highest-scoring, or longest window depending on user preference.
///
/// The implementation walks consecutive runs of identical
/// `(chrom,start,end,group)` windows and selects a representative. Score-based
/// policies fall back to the first instance when scores are absent.
///
/// # Parameters
/// - `windows`: sorted windows to deduplicate.
/// - `policy`: deduplication strategy.
/// - `use_score`: whether score data is available.
///
/// # Returns
/// Deduplicated window vector.
pub fn deduplicate_identical(
    windows: Vec<FinalWindow>,
    policy: DedupKeep,
    use_score: bool,
) -> Vec<FinalWindow> {
    if windows.is_empty() || matches!(policy, DedupKeep::None) {
        return windows;
    }

    let mut result: Vec<FinalWindow> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut j = i + 1;
        while j < windows.len()
            && windows[j].chrom == windows[i].chrom
            && windows[j].start == windows[i].start
            && windows[j].end == windows[i].end
            && windows[j].group == windows[i].group
        {
            j += 1;
        }

        if j == i + 1 {
            result.push(windows[i].clone());
        } else {
            let keep_index = match policy {
                DedupKeep::None => unreachable!(),
                DedupKeep::KeepFirst => i,
                DedupKeep::KeepHighestScore => {
                    if !use_score {
                        i
                    } else {
                        let mut best = i;
                        let mut best_score = windows[i].score.unwrap_or(f32::MIN);
                        for k in (i + 1)..j {
                            let sc = windows[k].score.unwrap_or(f32::MIN);
                            if sc > best_score {
                                best_score = sc;
                                best = k;
                            }
                        }
                        best
                    }
                }
                DedupKeep::KeepLowestScore => {
                    if !use_score {
                        i
                    } else {
                        let mut best = i;
                        let mut best_score = windows[i].score.unwrap_or(f32::MAX);
                        for k in (i + 1)..j {
                            let sc = windows[k].score.unwrap_or(f32::MAX);
                            if sc < best_score {
                                best_score = sc;
                                best = k;
                            }
                        }
                        best
                    }
                }
                DedupKeep::KeepLongest => {
                    let mut best = i;
                    let mut best_len = windows[i].length();
                    for k in (i + 1)..j {
                        let len = windows[k].length();
                        if len > best_len {
                            best_len = len;
                            best = k;
                        }
                    }
                    best
                }
            };

            result.push(windows[keep_index].clone());
        }

        i = j;
    }
    result
}

/// Enforce minimum distance between successive windows in a group.
///
/// Ensures windows within the same `(group, chrom)` pair appear at least
/// `min_distance_bp` apart, resolving conflicts with deterministic policies.
///
/// The function consumes already sorted input, collecting overlapping windows in
/// `pending_by_key`. Once spacing is satisfied, the chosen candidate is emitted
/// and the last end coordinate updated.
///
/// # Parameters
/// - `windows`: sorted windows to filter.
/// - `min_distance_bp`: required spacing (if any).
/// - `policy`: tie-breaking rule.
/// - `use_score`: indicates whether the score field is populated.
///
/// # Returns
/// Spaced window vector.
pub fn enforce_min_distance_within_group(
    windows: Vec<FinalWindow>,
    min_distance_bp: Option<u32>,
    policy: DistanceTiesPolicy,
    use_score: bool,
) -> Vec<FinalWindow> {
    if min_distance_bp.is_none() || windows.is_empty() {
        return windows;
    }
    let limit = min_distance_bp.unwrap();

    let mut result: Vec<FinalWindow> = Vec::with_capacity(windows.len());
    let mut last_end_by_key: FxHashMap<(String, String), u32> = FxHashMap::default();
    let mut pending_by_key: FxHashMap<(String, String), Vec<Candidate>> = FxHashMap::default();

    for (idx, window) in windows.iter().enumerate() {
        let key = (window.group.clone(), window.chrom.clone());
        let last_end = *last_end_by_key.get(&key).unwrap_or(&0);
        if window.start >= last_end.saturating_add(limit) {
            if let Some(cands) = pending_by_key.remove(&key) {
                let chosen = choose_candidate(&cands, policy, use_score);
                result.push(windows[chosen.window_idx].clone());
                last_end_by_key.insert(key.clone(), windows[chosen.window_idx].end);
            }
            pending_by_key.insert(
                key.clone(),
                vec![Candidate {
                    window_idx: idx,
                    score: window.score,
                    length: window.length(),
                }],
            );
        } else {
            pending_by_key
                .entry(key.clone())
                .or_default()
                .push(Candidate {
                    window_idx: idx,
                    score: window.score,
                    length: window.length(),
                });
        }
    }

    for (key, cands) in pending_by_key.into_iter() {
        let chosen = choose_candidate(&cands, policy, use_score);
        let w = windows[chosen.window_idx].clone();
        last_end_by_key.insert(key, w.end);
        result.push(w);
    }

    result
}

/// Split processed windows into a safe prefix and boundary tail.
///
/// Streaming emits the safe prefix immediately while retaining the tail to
/// protect against merges or spacing interactions with the next chunk.
///
/// The function computes the earliest index that might participate in a future
/// interaction by scanning backwards with the appropriate margin for the
/// current merge scope.
///
/// # Parameters
/// - `windows`: processed, sorted windows.
/// - `min_distance_bp`: optional spacing constraint.
/// - `merge_scope`: merge strategy in effect.
/// - `merge_gap_bp`: configured merge gap (if any).
///
/// # Returns
/// Tuple of `(safe_prefix, boundary_tail)`.
pub fn partition_safe_and_tail(
    windows: Vec<FinalWindow>,
    min_distance_bp: Option<u32>,
    merge_scope: MergeScope,
    merge_gap_bp: Option<u32>,
) -> (Vec<FinalWindow>, Vec<FinalWindow>) {
    if windows.is_empty() {
        return (windows, Vec::new());
    }
    let base_margin = min_distance_bp.unwrap_or(0).max(merge_gap_bp.unwrap_or(0));

    let mut tail_start_index = windows.len();
    let mut expanded = true;

    while expanded {
        expanded = false;
        let candidate_index = match merge_scope {
            MergeScope::None => windows.len(),
            MergeScope::Within => compute_tail_start_within(&windows, base_margin),
            MergeScope::Across => compute_tail_start_across(&windows, base_margin),
        };

        if candidate_index < tail_start_index {
            tail_start_index = candidate_index;
            expanded = true;
        }
    }

    let safe_prefix = windows[..tail_start_index].to_vec();
    let boundary_tail = windows[tail_start_index..].to_vec();
    (safe_prefix, boundary_tail)
}

/// Identify earliest index in the suffix that might be affected within groups.
///
/// Walks backwards tracking the most recent window per `(group, chrom)` to see
/// how far into the suffix future merges or spacing constraints could reach.
///
/// Keeps a map of `last_end` per key and expands the unsafe region whenever a
/// window overlaps the allowed margin.
///
/// # Parameters
/// - `windows`: processed windows under consideration.
/// - `margin`: distance margin derived from spacing/merge settings.
///
/// # Returns
/// Earliest index that must remain in the tail.
fn compute_tail_start_within(windows: &[FinalWindow], margin: u32) -> usize {
    let mut last_per_key: FxHashMap<(String, String), u32> = FxHashMap::default();
    let mut min_index = windows.len();

    for (idx, window) in windows.iter().enumerate().rev() {
        let key = (window.group.clone(), window.chrom.clone());
        let maybe_last_end = last_per_key.get(&key).copied();
        match maybe_last_end {
            None => {
                last_per_key.insert(key, window.end);
                if margin > 0 {
                    min_index = min_index.min(idx);
                }
            }
            Some(last_end) => {
                if window.start <= last_end.saturating_add(margin) {
                    min_index = min_index.min(idx);
                    let new_end = last_end.max(window.end);
                    last_per_key.insert(key, new_end);
                }
            }
        }
    }
    min_index
}

/// Identify earliest index in the suffix that might be affected across groups.
///
/// Similar to `compute_tail_start_within` but tracked purely by chromosome when
/// merging across groups.
///
/// Maintains the last end coordinate per chromosome and expands the unsafe
/// region when the margin is violated.
///
/// # Parameters
/// - `windows`: processed windows under consideration.
/// - `margin`: distance margin derived from spacing/merge settings.
///
/// # Returns
/// Earliest index that must remain in the tail.
fn compute_tail_start_across(windows: &[FinalWindow], margin: u32) -> usize {
    let mut last_end_by_chrom: FxHashMap<&str, u32> = FxHashMap::default();
    let mut min_index = windows.len();

    for (idx, window) in windows.iter().enumerate().rev() {
        let last_end = last_end_by_chrom
            .get(window.chrom.as_str())
            .copied()
            .unwrap_or(0);
        if last_end == 0 {
            last_end_by_chrom.insert(window.chrom.as_str(), window.end);
            if margin > 0 {
                min_index = min_index.min(idx);
            }
        } else if window.start <= last_end.saturating_add(margin) {
            min_index = min_index.min(idx);
            let new_end = last_end.max(window.end);
            last_end_by_chrom.insert(window.chrom.as_str(), new_end);
        }
    }
    min_index
}

fn choose_candidate(
    candidates: &[Candidate],
    policy: DistanceTiesPolicy,
    use_score: bool,
) -> Candidate {
    match policy {
        DistanceTiesPolicy::KeepFirst => candidates[0].clone(),
        DistanceTiesPolicy::KeepHighestScore => {
            if use_score {
                candidates
                    .iter()
                    .max_by(|a, b| {
                        a.score
                            .unwrap_or(f32::MIN)
                            .partial_cmp(&b.score.unwrap_or(f32::MIN))
                            .unwrap_or(Ordering::Equal)
                    })
                    .unwrap()
                    .clone()
            } else {
                candidates[0].clone()
            }
        }
        DistanceTiesPolicy::KeepLowestScore => {
            if use_score {
                candidates
                    .iter()
                    .min_by(|a, b| {
                        a.score
                            .unwrap_or(f32::MAX)
                            .partial_cmp(&b.score.unwrap_or(f32::MAX))
                            .unwrap_or(Ordering::Equal)
                    })
                    .unwrap()
                    .clone()
            } else {
                candidates[0].clone()
            }
        }
        DistanceTiesPolicy::KeepLongest => {
            candidates.iter().max_by_key(|c| c.length).unwrap().clone()
        }
    }
}
