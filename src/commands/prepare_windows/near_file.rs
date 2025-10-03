use crate::shared::io::open_text_reader;
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use std::{cmp::Ordering, io::BufRead, path::Path};

use crate::commands::prepare_windows::config::NearEdge;

/// Interval from the `--near` set.
///
/// The `group` is optional; it is used when composing output groups with distance bins.
#[derive(Debug, Clone, Copy)]
pub struct NearInterval {
    pub start: u32,
    pub end: u32,
    pub group_id: Option<u32>,
}

/// Per-chromosome near intervals, sorted and validated.
#[derive(Debug, Default)]
pub struct NearChrom {
    pub intervals: Vec<NearInterval>,
}

/// Index holding near intervals per chromosome and a compact group-id interner.
#[derive(Debug, Default)]
pub struct NearIndex {
    pub per_chrom: FxHashMap<String, NearChrom>,
    pub group_name_to_id: FxHashMap<String, u32>,
    pub group_id_to_name: Vec<String>,
}

/// Where the nearest interval sits relative to the window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NearSide {
    Upstream,
    Downstream,
    Overlap,
}

/// Distance, group id, and side for a single nearest hit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NearHit {
    pub distance: i32,
    pub group_id: Option<u32>,
    pub side: NearSide,
}

/// Captures the upstream and downstream hits when a tie happens.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NearTie {
    pub upstream: Option<NearHit>,
    pub downstream: Option<NearHit>,
}

/// Result of the nearest-edge lookup: either one hit or a tie.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NearestResult {
    Single(NearHit),
    Tie(NearTie),
}

/// Load and validate the `--near` intervals.
///
/// Validation rules:
/// - Intervals must have end > start.
/// - Intervals must be pairwise non-overlapping per chromosome (half-open).
/// - No two intervals on the same chromosome may have identical (start, end) edges.
/// These rules avoid ambiguous "nearest" selections and keep distance semantics clear.
///
/// Parameters
/// ----------
/// - path:
///     Path to the near file.
/// - separator:
///     Field separator character.
/// - has_header:
///     Whether the near file has a header line.
/// - group_col_present:
///     Whether a `group` column is present (fourth column).
///
/// Returns
/// -------
/// - index:
///     Per-chromosome near intervals and a group interner.
pub fn load_near_index(
    path: &Path,
    separator: char,
    has_header: bool,
    group_col_present: bool,
) -> Result<NearIndex> {
    let mut reader = open_text_reader(path)?;
    let mut line = String::new();

    if has_header {
        line.clear();
        reader.read_line(&mut line)?; // discard header
    }

    let mut raw_by_chrom: FxHashMap<String, Vec<(u32, u32, Option<String>)>> = FxHashMap::default();

    for (lineno, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split(separator);

        let chrom = it
            .next()
            .context("Near parse: missing chrom")?
            .trim()
            .to_string();
        let start: u32 = it
            .next()
            .context("Near parse: missing start")?
            .trim()
            .parse()
            .with_context(|| format!("Invalid start at near line {}", lineno + 1))?;
        let end: u32 = it
            .next()
            .context("Near parse: missing end")?
            .trim()
            .parse()
            .with_context(|| format!("Invalid end at near line {}", lineno + 1))?;

        if end <= start {
            bail!(
                "Near parse error at line {}: end ({}) must be > start ({})",
                lineno + 1,
                end,
                start
            );
        }

        let group = if group_col_present {
            it.next()
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
        } else {
            None
        };

        raw_by_chrom
            .entry(chrom)
            .or_default()
            .push((start, end, group));
    }

    // Validate and intern groups
    let mut index = NearIndex::default();

    for (chrom, mut items) in raw_by_chrom {
        // Sort by (start, end, group)
        items.sort_unstable_by(|a, b| {
            a.0.cmp(&b.0)
                .then(a.1.cmp(&b.1))
                .then_with(|| match (&a.2, &b.2) {
                    (Some(ga), Some(gb)) => ga.cmp(gb),
                    (None, Some(_)) => Ordering::Less,
                    (Some(_), None) => Ordering::Greater,
                    (None, None) => Ordering::Equal,
                })
        });

        // Validate no overlap and no identical edges per chromosome
        let mut validated: Vec<NearInterval> = Vec::with_capacity(items.len());
        let mut last_end: u32 = 0;
        let mut last_start: u32 = 0;
        let mut have_last = false;

        for (start, end, group_opt) in items.into_iter() {
            if have_last {
                if start < last_end {
                    bail!(
                        "Near validation failed on {}: intervals overlap at [{}, {}) and previous ending at {}",
                        chrom,
                        start,
                        end,
                        last_end
                    );
                }
                if start == last_start && end == last_end {
                    bail!(
                        "Near validation failed on {}: duplicate edges at [{}, {})",
                        chrom,
                        start,
                        end
                    );
                }
            }
            last_start = start;
            last_end = end;
            have_last = true;

            let group_id = if let Some(name) = group_opt {
                Some(
                    *index
                        .group_name_to_id
                        .entry(name.clone())
                        .or_insert_with(|| {
                            let id = index.group_id_to_name.len() as u32;
                            index.group_id_to_name.push(name);
                            id
                        }),
                )
            } else {
                None
            };

            validated.push(NearInterval {
                start,
                end,
                group_id,
            });
        }

        index.per_chrom.insert(
            chrom,
            NearChrom {
                intervals: validated,
            },
        );
    }

    Ok(index)
}

/// Compute distance from a window to the nearest `near` interval edge.
///
/// Concept:
/// - If the window overlaps the interval, the distance is 0 (we are at the site).
/// - Else, the distance is the minimum of the distances from the window's two edges
///   to the selected target edges (left, right, or nearest).
///
/// Signed distances:
/// - Upstream is negative, downstream is positive, relative to the chosen nearest interval.
///
/// Parameters
/// ----------
/// - window_start:
///     Window start.
/// - window_end:
///     Window end (exclusive).
/// - near_chrom:
///     Validated, sorted intervals on the same chromosome.
/// - which_edge:
///     Edge selection mode (Left, Right, Nearest).
/// - signed:
///     Whether to return signed distances.
///
/// Returns
/// -------
/// - result:
///     Either a single nearest hit or a tie describing upstream/downstream hits.
pub fn nearest_edge_distance(
    window_start: u32,
    window_end: u32,
    near_chrom: &NearChrom,
    which_edge: &NearEdge,
    signed: bool,
) -> Option<NearestResult> {
    if near_chrom.intervals.is_empty() {
        return None;
    }

    // Track the best absolute distance and every interval that matches it.
    #[derive(Clone, Copy)]
    struct Candidate {
        idx: usize,
        distance: i32,
        side: NearSide,
    }

    let idx = near_chrom
        .intervals
        .binary_search_by_key(&window_start, |iv| iv.start)
        .unwrap_or_else(|i| i);

    let mut best_abs_distance: Option<i32> = None;
    let mut best_candidates: Vec<Candidate> = Vec::new();

    let start_index = idx.saturating_sub(2);
    let end_index = (idx + 2).min(near_chrom.intervals.len().saturating_sub(1));

    for j in start_index..=end_index {
        let iv = near_chrom.intervals[j];

        if window_start < iv.end && window_end > iv.start {
            return Some(NearestResult::Single(NearHit {
                distance: 0,
                group_id: iv.group_id,
                side: NearSide::Overlap,
            }));
        }

        let mut candidate_edges: [i32; 2] = [iv.start as i32, iv.end as i32];
        let num_edges = match which_edge {
            NearEdge::Left => 1,
            NearEdge::Right => {
                candidate_edges[0] = iv.end as i32;
                1
            }
            NearEdge::Nearest => 2,
        };

        for k in 0..num_edges {
            let target_edge = candidate_edges[k];
            let diff_start = target_edge - window_start as i32;
            let diff_end = target_edge - window_end as i32;

            let mut signed_distance = if diff_start.abs() <= diff_end.abs() {
                diff_start
            } else {
                diff_end
            };

            // Direction follows the sign of the chosen difference: negative means
            // the interval lies upstream of the window, positive means downstream.
            let side = if signed_distance < 0 {
                NearSide::Upstream
            } else if signed_distance > 0 {
                NearSide::Downstream
            } else {
                NearSide::Overlap
            };

            let abs_distance = signed_distance.abs();
            if !signed {
                signed_distance = abs_distance;
            }

            match best_abs_distance {
                None => {
                    best_abs_distance = Some(abs_distance);
                    best_candidates.clear();
                    best_candidates.push(Candidate {
                        idx: j,
                        distance: signed_distance,
                        side,
                    });
                }
                Some(current_abs) => {
                    if abs_distance < current_abs {
                        best_abs_distance = Some(abs_distance);
                        best_candidates.clear();
                        best_candidates.push(Candidate {
                            idx: j,
                            distance: signed_distance,
                            side,
                        });
                    } else if abs_distance == current_abs {
                        if !best_candidates.iter().any(|cand| cand.idx == j) {
                            best_candidates.push(Candidate {
                                idx: j,
                                distance: signed_distance,
                                side,
                            });
                        }
                    }
                }
            }
        }
    }

    if best_candidates.is_empty() {
        return None;
    }

    if best_candidates.len() == 1 {
        let cand = best_candidates[0];
        let group_id = near_chrom.intervals[cand.idx].group_id;
        return Some(NearestResult::Single(NearHit {
            distance: cand.distance,
            group_id,
            side: cand.side,
        }));
    }

    // More than one interval shares the same distance: record the upstream and
    // downstream hits so the caller can decide how to label or drop them.
    let mut tie = NearTie::default();
    for cand in best_candidates {
        let group_id = near_chrom.intervals[cand.idx].group_id;
        let hit = NearHit {
            distance: cand.distance,
            group_id,
            side: cand.side,
        };
        match cand.side {
            NearSide::Upstream => {
                if tie.upstream.is_none() {
                    tie.upstream = Some(hit);
                }
            }
            NearSide::Downstream => {
                if tie.downstream.is_none() {
                    tie.downstream = Some(hit);
                }
            }
            NearSide::Overlap => {
                return Some(NearestResult::Single(hit));
            }
        }
    }

    Some(NearestResult::Tie(tie))
}
