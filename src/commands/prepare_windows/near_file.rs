use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use std::{
    cmp::Ordering,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

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
    let file = File::open(path).with_context(|| format!("Opening near file {:?}", path))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
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
/// - distance:
///     The distance (0 if overlapping).
/// - near_group_id:
///     The near interval's group id if a unique nearest interval exists; otherwise None.
pub fn nearest_edge_distance(
    window_start: u32,
    window_end: u32,
    near_chrom: &NearChrom,
    which_edge: &NearEdge,
    signed: bool,
) -> Option<(i32, Option<u32>)> {
    if near_chrom.intervals.is_empty() {
        return None;
    }

    // Binary search candidate by start
    let idx = near_chrom
        .intervals
        .binary_search_by_key(&window_start, |iv| iv.start)
        .unwrap_or_else(|i| i);

    // We will check a small neighborhood around idx. Because the near intervals are
    // non-overlapping and sorted, the nearest interval (by edge distance) must be among
    // a small set near the insertion point.
    let mut best_distance: Option<i32> = None;
    let mut best_group_id: Option<u32> = None;
    let mut best_interval_index: Option<usize> = None;

    let start_index = idx.saturating_sub(2);
    let end_index = (idx + 2).min(near_chrom.intervals.len().saturating_sub(1));

    for j in start_index..=end_index {
        let iv = near_chrom.intervals[j];

        // Overlap check: window overlaps interval -> distance = 0
        if window_start < iv.end && window_end > iv.start {
            return Some((0, iv.group_id));
        }

        // Candidate target edges
        let mut candidate_edges: [i32; 2] = [iv.start as i32, iv.end as i32]; // right edge is at end (exclusive)
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

            // Distance to this edge from window edges
            let distance_from_start = (window_start as i32) - target_edge;
            let distance_from_end = (window_end as i32) - target_edge;

            // Choose the closer window edge for this target edge
            let mut distance = if distance_from_start.abs() <= distance_from_end.abs() {
                distance_from_start
            } else {
                distance_from_end
            };

            if !signed {
                distance = distance.abs();
            }

            match best_distance {
                None => {
                    best_distance = Some(distance);
                    best_group_id = iv.group_id;
                    best_interval_index = Some(j);
                }
                Some(current) => {
                    if distance.abs() < current.abs() {
                        best_distance = Some(distance);
                        best_group_id = iv.group_id;
                        best_interval_index = Some(j);
                    } else if distance.abs() == current.abs() {
                        // Tie between different intervals or edges. Because we validated that
                        // intervals do not overlap and do not share identical edges, an exact
                        // tie can only occur in gaps symmetrical to two neighbors.
                        //
                        // Design choice to avoid arbitrary selection:
                        // - Keep the distance value as-is (absolute tie).
                        // - Do not attach a near group label in this tie case (set to None).
                        // - If signed is requested, sign will follow the interval with smaller start
                        //   deterministically below.
                        best_group_id = None;

                        // Resolve sign deterministically if needed by picking the left interval.
                        if signed {
                            if let Some(prev_idx) = best_interval_index {
                                let left_idx = j.min(prev_idx);
                                let left_iv = near_chrom.intervals[left_idx];
                                // Recompute sign vs left interval's nearest edge
                                let left_edge = if (window_start as i32 - left_iv.end as i32).abs()
                                    <= (window_end as i32 - left_iv.end as i32).abs()
                                {
                                    left_iv.end as i32
                                } else {
                                    left_iv.start as i32
                                };
                                let signed_from_start = (window_start as i32) - left_edge;
                                let signed_from_end = (window_end as i32) - left_edge;
                                let signed_distance =
                                    if signed_from_start.abs() <= signed_from_end.abs() {
                                        signed_from_start
                                    } else {
                                        signed_from_end
                                    };
                                best_distance = Some(signed_distance);
                                best_interval_index = Some(left_idx);
                            }
                        } else {
                            best_distance = Some(distance.abs());
                        }
                    }
                }
            }
        }
    }

    best_distance.map(|d| (d, best_group_id))
}
