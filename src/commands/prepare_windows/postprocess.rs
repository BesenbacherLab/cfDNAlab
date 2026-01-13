use crate::commands::prepare_windows::{
    config::{CoordinateSet, DedupKeep, DistancePolicy, MergeScope},
    prepare_windows::Window,
};
use fxhash::FxHashMap;
use std::cmp::Ordering;
use std::sync::Arc;

/// Candidate retained while resolving minimum-distance conflicts.
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
/// `(chrom,start,end,group_key)` windows and selects a representative. Score-based
/// policies fall back to the first instance when scores are absent.
///
/// # Parameters
/// - `windows`: sorted windows to deduplicate.
/// - `policy`: deduplication strategy.
/// - `use_score`: whether score data is available.
/// - `coord_set`: coordinate set used to detect duplicates and compare lengths.
///
/// # Returns
/// Deduplicated window vector.
pub fn deduplicate_identical(
    windows: Vec<Window>,
    policy: DedupKeep,
    use_score: bool,
    coord_set: CoordinateSet,
) -> Vec<Window> {
    if windows.is_empty() || matches!(policy, DedupKeep::None) {
        return windows;
    }

    let mut result: Vec<Window> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut j = i + 1;
        while j < windows.len()
            && windows[j].chrom == windows[i].chrom
            && windows[j].start_for(coord_set) == windows[i].start_for(coord_set)
            && windows[j].end_for(coord_set) == windows[i].end_for(coord_set)
            && windows[j].group_key == windows[i].group_key
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
            };

            result.push(windows[keep_index].clone());
        }

        i = j;
    }
    result
}

/// Enforce minimum distance between successive windows in a group.
///
/// Ensures windows within the same `(group_key, chrom)` pair appear at least
/// `min_distance_bp` apart, resolving conflicts with deterministic policies.
///
/// The function consumes already sorted input, collecting overlapping windows in
/// `pending_by_key`. Once the minimum distance is satisfied, it selects the best candidate,
/// appends it to the result, and updates the last end coordinate.
///
/// # Parameters
/// - `windows`: sorted windows to filter.
/// - `min_distance_bp`: required minimum distance, or None to disable minimum-distance filtering.
/// - `policy`: tie-breaking rule.
/// - `use_score`: indicates whether the score field is populated.
/// - `coord_set`: coordinate set used to compute distances and lengths.
///
/// # Returns
/// Spaced window vector.
pub fn enforce_min_distance_within_group(
    windows: Vec<Window>,
    min_distance_bp: Option<u32>,
    policy: DistancePolicy,
    use_score: bool,
    coord_set: CoordinateSet,
) -> Vec<Window> {
    if min_distance_bp.is_none() || windows.is_empty() {
        return windows;
    }
    let limit = min_distance_bp.unwrap();

    let mut result: Vec<Window> = Vec::with_capacity(windows.len());
    let mut last_end_by_key: FxHashMap<(String, Arc<str>), u32> = FxHashMap::default();
    let mut pending_by_key: FxHashMap<(String, Arc<str>), Vec<Candidate>> = FxHashMap::default();

    for (idx, window) in windows.iter().enumerate() {
        let key = (window.group_key.clone(), window.chrom.clone());
        let last_end = *last_end_by_key.get(&key).unwrap_or(&0);
        if window.start_for(coord_set) >= last_end.saturating_add(limit) {
            if let Some(cands) = pending_by_key.remove(&key) {
                let chosen = choose_candidate(&cands, policy, use_score);
                result.push(windows[chosen.window_idx].clone());
                last_end_by_key.insert(key.clone(), windows[chosen.window_idx].end_for(coord_set));
            }
            pending_by_key.insert(
                key.clone(),
                vec![Candidate {
                    window_idx: idx,
                    score: window.score,
                    length: window.length_for(coord_set),
                }],
            );
        } else {
            pending_by_key
                .entry(key.clone())
                .or_default()
                .push(Candidate {
                    window_idx: idx,
                    score: window.score,
                    length: window.length_for(coord_set),
                });
        }
    }

    for (key, cands) in pending_by_key.into_iter() {
        let chosen = choose_candidate(&cands, policy, use_score);
        let w = windows[chosen.window_idx].clone();
        last_end_by_key.insert(key, w.end_for(coord_set));
        result.push(w);
    }

    result
}

/// Tag windows as clustered based on average position-wise overlap across groups.
///
/// This computes the average overlap depth within each chromosome and assigns
/// `cluster=cluster` when it meets the threshold, otherwise `cluster=none`.
/// The average is the total overlap depth across the window divided by its length.
/// The input must already be sorted by `(chrom, start, end)` in the coordinate
/// set used for overlap checks.
///
/// Parameters
/// ----------
/// - `windows`:
///     Sorted windows to label in place.
///
/// - `min_overlaps`:
///     Minimum average overlap depth required, counting the window itself.
///
/// - `coord_set`:
///     Coordinate set used for overlap checks.
///
/// Returns
/// -------
/// - `()`:
///     Updates `windows` in place.
pub fn apply_cluster_labels(windows: &mut [Window], min_overlaps: u32, coord_set: CoordinateSet) {
    if windows.is_empty() {
        return;
    }

    let min_avg_overlap = min_overlaps as u64;
    let cluster_label = "cluster".to_string();
    let non_cluster_label = "none".to_string();

    let mut start_idx = 0usize;
    while start_idx < windows.len() {
        let chrom = windows[start_idx].chrom.clone();
        let mut end_idx = start_idx + 1;
        while end_idx < windows.len() && windows[end_idx].chrom == chrom {
            end_idx += 1;
        }

        let (boundaries, coverage_by_segment, coverage_prefix) =
            build_coverage_index(&windows[start_idx..end_idx], coord_set);

        if boundaries.len() < 2 {
            for window in &mut windows[start_idx..end_idx] {
                for tuple in &mut window.label_tuples {
                    tuple.cluster = Some(non_cluster_label.clone());
                }
            }
            start_idx = end_idx;
            continue;
        }

        let window_count = end_idx - start_idx;
        let mut start_segment_by_window = vec![0usize; window_count];

        // Map each window start to the segment containing it using a linear sweep
        let mut segment_idx = 0usize;
        for (offset, window) in windows[start_idx..end_idx].iter().enumerate() {
            let window_start = window.start_for(coord_set);
            while segment_idx + 1 < boundaries.len() && boundaries[segment_idx + 1] <= window_start
            {
                segment_idx += 1;
            }
            start_segment_by_window[offset] = segment_idx;
        }

        for (offset, window) in windows[start_idx..end_idx].iter_mut().enumerate() {
            let window_start = window.start_for(coord_set);
            let window_end = window.end_for(coord_set);
            let window_length = window_end.saturating_sub(window_start) as u64;
            // Find the end segment by stepping forward from the start segment
            let mut end_segment_idx = start_segment_by_window[offset];
            while end_segment_idx + 1 < boundaries.len()
                && boundaries[end_segment_idx + 1] < window_end
            {
                end_segment_idx += 1;
            }
            // Compute total overlap depth by summing coverage across the window span
            let total_overlap = overlap_sum_with_segments(
                window_start,
                window_end,
                start_segment_by_window[offset],
                end_segment_idx,
                &boundaries,
                &coverage_by_segment,
                &coverage_prefix,
            );
            let is_cluster = total_overlap >= min_avg_overlap.saturating_mul(window_length);
            let label = if is_cluster {
                &cluster_label
            } else {
                &non_cluster_label
            };

            for tuple in &mut window.label_tuples {
                tuple.cluster = Some(label.clone());
            }
        }

        start_idx = end_idx;
    }
}

fn build_coverage_index(
    windows: &[Window],
    coord_set: CoordinateSet,
) -> (Vec<u32>, Vec<u32>, Vec<u64>) {
    let mut events: Vec<(u32, i32)> = Vec::with_capacity(windows.len() * 2);
    for window in windows {
        events.push((window.start_for(coord_set), 1));
        events.push((window.end_for(coord_set), -1));
    }
    events.sort_unstable_by(|a, b| a.0.cmp(&b.0));

    let mut boundaries: Vec<u32> = Vec::new();
    let mut deltas: Vec<i32> = Vec::new();
    for (pos, delta) in events {
        if boundaries.last().copied() == Some(pos) {
            let slot = deltas.last_mut().expect("delta slot");
            *slot += delta;
        } else {
            boundaries.push(pos);
            deltas.push(delta);
        }
    }

    let mut coverage_by_segment: Vec<u32> = Vec::with_capacity(boundaries.len().saturating_sub(1));
    let mut coverage_prefix: Vec<u64> = Vec::with_capacity(boundaries.len());
    coverage_prefix.push(0);

    let mut coverage_depth: i32 = 0;
    for i in 0..boundaries.len().saturating_sub(1) {
        coverage_depth += deltas[i];
        let segment_len = boundaries[i + 1].saturating_sub(boundaries[i]);
        let segment_coverage = coverage_depth.max(0) as u32;
        coverage_by_segment.push(segment_coverage);
        let cumulative = coverage_prefix[i] + (segment_coverage as u64) * (segment_len as u64);
        coverage_prefix.push(cumulative);
    }

    (boundaries, coverage_by_segment, coverage_prefix)
}

fn overlap_sum_with_segments(
    start: u32,
    end: u32,
    start_segment_idx: usize,
    end_segment_idx: usize,
    boundaries: &[u32],
    coverage_by_segment: &[u32],
    coverage_prefix: &[u64],
) -> u64 {
    if start >= end || boundaries.len() < 2 {
        return 0;
    }

    let max_segment = coverage_by_segment.len().saturating_sub(1);
    let start_idx = start_segment_idx.min(max_segment);
    let end_idx = end_segment_idx.min(max_segment);

    if start_idx == end_idx {
        let segment_len = end.saturating_sub(start) as u64;
        return (coverage_by_segment[start_idx] as u64) * segment_len;
    }

    let mut total = 0u64;
    let first_len = boundaries[start_idx + 1].saturating_sub(start) as u64;
    total += (coverage_by_segment[start_idx] as u64) * first_len;

    if end_idx > start_idx + 1 {
        total += coverage_prefix[end_idx] - coverage_prefix[start_idx + 1];
    }

    let last_len = end.saturating_sub(boundaries[end_idx]) as u64;
    total += (coverage_by_segment[end_idx] as u64) * last_len;
    total
}

/// Split processed windows into a processed region and boundary tail.
///
/// The processed region cannot be affected by any future windows, so it can be
/// written immediately. The tail stays in memory in case future windows
/// merge, overlap, or fall within the minimum-distance window.
///
/// The function finds the earliest output-order index that could still interact
/// across the chunk boundary under either merging or minimum-distance rules, then splits
/// the list.
/// When a margin is zero, the last window for each key or chromosome is still kept in the tail
/// because overlap is still possible across the chunk boundary.
///
/// # Parameters
/// - `windows`: processed windows in output order.
/// - `min_distance_bp`: optional minimum-distance constraint.
/// - `merge_scope`: merge strategy in effect.
/// - `merge_gap_bp`: configured merge gap (if any).
/// - `merge_coord_set`: coordinate set used for merging.
/// - `distance_coord_set`: coordinate set used for minimum-distance checks.
/// - `cluster_min_overlaps`: optional overlap threshold for cluster labeling.
/// - `cluster_coord_set`: coordinate set used for cluster overlap checks.
/// - `merge_group_keys`: optional merge-group values aligned with `windows`.
///
/// # Returns
/// Tuple of `(safe_prefix, boundary_tail)`.
pub fn partition_safe_and_tail(
    windows: Vec<Window>,
    min_distance_bp: Option<u32>,
    merge_scope: MergeScope,
    merge_gap_bp: Option<u32>,
    merge_coord_set: CoordinateSet,
    distance_coord_set: CoordinateSet,
    cluster_min_overlaps: Option<u32>,
    cluster_coord_set: CoordinateSet,
    merge_group_keys: Option<&[String]>,
) -> (Vec<Window>, Vec<Window>) {
    if windows.is_empty() {
        return (windows, Vec::new());
    }

    // Index of the first window that might still interact with a later chunk
    let mut tail_start_index = windows.len();

    // Merging can pull later windows leftward into the current chunk
    if let Some(gap_bp) = merge_gap_bp {
        if !matches!(merge_scope, MergeScope::None) {
            if let Some(min_index) = collect_tail_indices(
                &windows,
                merge_scope,
                gap_bp,
                merge_coord_set,
                merge_group_keys,
            )
            .into_iter()
            .min()
            {
                tail_start_index = tail_start_index.min(min_index);
                println!("merge tail start index: {}", min_index);
            }
        }
    }

    // Minimum-distance filtering can also skip across the chunk boundary
    if let Some(distance_bp) = min_distance_bp {
        if let Some(min_index) = collect_tail_indices(
            &windows,
            MergeScope::Within,
            distance_bp,
            distance_coord_set,
            None,
        )
        .into_iter()
        .min()
        {
            tail_start_index = tail_start_index.min(min_index);
            println!("min-distance tail start index: {}", min_index);
        }
    }

    // Clustering depends on overlap depth, so keep any window that could overlap
    if let Some(min_overlaps) = cluster_min_overlaps {
        if min_overlaps > 1 {
            // For overlap-only clustering (margin is zero), any window whose end is <= the first
            // start of the yet-to-be-seen chunk cannot overlap future windows. Everything that
            // extends beyond that boundary must stay in the tail.
            if let Some(max_start) = windows.iter().map(|w| w.start_for(cluster_coord_set)).max() {
                if let Some(min_index) = windows.iter().enumerate().find_map(|(idx, w)| {
                    if w.end_for(cluster_coord_set) > max_start {
                        Some(idx)
                    } else {
                        None
                    }
                }) {
                    tail_start_index = tail_start_index.min(min_index);
                    println!("clustering tail start index: {}", min_index);
                }
            }
        }
    }

    // If no risk was found, everything is safe to write
    if tail_start_index == windows.len() {
        return (windows, Vec::new());
    }

    let safe_prefix = windows[..tail_start_index].to_vec();
    let boundary_tail = windows[tail_start_index..].to_vec();
    (safe_prefix, boundary_tail)
}

fn collect_tail_indices(
    windows: &[Window],
    merge_scope: MergeScope,
    margin: u32,
    coord_set: CoordinateSet,
    group_keys: Option<&[String]>,
) -> Vec<usize> {
    if windows.is_empty() {
        return Vec::new();
    }

    let mut indices: Vec<usize> = (0..windows.len()).collect();
    match merge_scope {
        MergeScope::Within => indices.sort_unstable_by(|&a_idx, &b_idx| {
            let a = &windows[a_idx];
            let b = &windows[b_idx];
            let a_key = group_keys
                .and_then(|keys| keys.get(a_idx))
                .map(|key| key.as_str())
                .unwrap_or_else(|| a.group_key.as_str());
            let b_key = group_keys
                .and_then(|keys| keys.get(b_idx))
                .map(|key| key.as_str())
                .unwrap_or_else(|| b.group_key.as_str());
            a_key
                .cmp(b_key)
                .then(a.chrom.cmp(&b.chrom))
                .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
        }),
        MergeScope::Across => indices.sort_unstable_by(|&a_idx, &b_idx| {
            let a = &windows[a_idx];
            let b = &windows[b_idx];
            a.chrom
                .cmp(&b.chrom)
                .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
                .then(a.group_key.cmp(&b.group_key))
        }),
        MergeScope::None => return Vec::new(),
    }

    let tail_start = match merge_scope {
        MergeScope::None => indices.len(),
        MergeScope::Within => {
            compute_tail_start_within_indices(&indices, windows, margin, coord_set, group_keys)
        }
        MergeScope::Across => {
            compute_tail_start_across_indices(&indices, windows, margin, coord_set)
        }
    };

    indices[tail_start..].to_vec()
}

/// Identify earliest index in the suffix that might be affected within groups.
///
/// Walks backwards tracking the most recent window per `(group_key, chrom)` to see
/// how far into the suffix future merges or minimum-distance filtering could reach.
///
/// Keeps a map of `last_end` per key and expands the unsafe region whenever a
/// window overlaps the allowed margin.
///
/// # Parameters
/// - `indices`: sorted window indices under consideration.
/// - `windows`: window storage for index lookups.
/// - `margin`: distance margin derived from minimum-distance and merge settings.
/// - `coord_set`: coordinate set used for overlap checks.
///
/// # Returns
/// Earliest index that must remain in the tail.
fn compute_tail_start_within_indices(
    indices: &[usize],
    windows: &[Window],
    margin: u32,
    coord_set: CoordinateSet,
    group_keys: Option<&[String]>,
) -> usize {
    // We walk the sorted indices from the end toward the start so we can see which earlier
    // windows could still be impacted by a later window in this chunk
    // For minimum-distance and within-group merges, a window is unsafe to flush if there exists
    // a later window in the same (group, chrom) whose start lies within `margin` of this window's start
    // We track the farthest end we have seen so far per key to decide whether the current window
    // overlaps or sits within the margin of that span
    let mut last_per_key: FxHashMap<(String, Arc<str>), u32> = FxHashMap::default();
    let mut min_index = indices.len();

    for (pos, &idx) in indices.iter().enumerate().rev() {
        // Current window when scanning backward
        let window = &windows[idx];
        // Resolve group key used for within-group rules
        let group_key = group_keys
            .and_then(|keys| keys.get(idx))
            .cloned()
            .unwrap_or_else(|| window.group_key.clone());
        let key = (group_key, window.chrom.clone());
        let maybe_last_end = last_per_key.get(&key).copied();
        let window_end = window.end_for(coord_set);
        match maybe_last_end {
            None => {
                // First time seeing this (group, chrom) when scanning from the end
                // Record its end so earlier windows can measure overlap against it
                last_per_key.insert(key, window_end);
            }
            Some(last_end) => {
                // There is already a later window for this key
                // If the current start is within the margin of that later end, this window is unsafe
                // and everything at or before this position must be kept in the tail
                if window.start_for(coord_set) <= last_end.saturating_add(margin) {
                    min_index = min_index.min(pos);
                    let new_end = last_end.max(window_end);
                    // Extend the tracked end so even earlier windows consider the union span
                    last_per_key.insert(key, new_end);
                }
            }
        }
    }
    min_index
}

/// Identify earliest index in the suffix that might be affected across groups.
///
/// Similar to `compute_tail_start_within_indices` but tracked purely by chromosome
/// when merging across groups.
///
/// Maintains the last end coordinate per chromosome and expands the unsafe
/// region when the margin is violated.
///
/// # Parameters
/// - `indices`: sorted window indices under consideration.
/// - `windows`: window storage for index lookups.
/// - `margin`: distance margin derived from minimum-distance and merge settings.
/// - `coord_set`: coordinate set used for overlap checks.
///
/// # Returns
/// Earliest index that must remain in the tail.
fn compute_tail_start_across_indices(
    indices: &[usize],
    windows: &[Window],
    margin: u32,
    coord_set: CoordinateSet,
) -> usize {
    let mut last_end_by_chrom: FxHashMap<&str, u32> = FxHashMap::default();
    let mut min_index = indices.len();

    for (pos, &idx) in indices.iter().enumerate().rev() {
        let window = &windows[idx];
        let chrom_name = window.chrom.as_ref();
        let last_end = last_end_by_chrom.get(chrom_name).copied().unwrap_or(0);
        let window_end = window.end_for(coord_set);

        if last_end == 0 {
            last_end_by_chrom.insert(chrom_name, window_end);
        } else if window.start_for(coord_set) <= last_end.saturating_add(margin) {
            min_index = min_index.min(pos);
            let new_end = last_end.max(window_end);
            last_end_by_chrom.insert(chrom_name, new_end);
        }
    }
    min_index
}

fn choose_candidate(
    candidates: &[Candidate],
    policy: DistancePolicy,
    use_score: bool,
) -> Candidate {
    match policy {
        DistancePolicy::KeepFirst => candidates[0].clone(),
        DistancePolicy::KeepHighestScore => {
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
        DistancePolicy::KeepLowestScore => {
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
        DistancePolicy::KeepLongest => candidates.iter().max_by_key(|c| c.length).unwrap().clone(),
    }
}
