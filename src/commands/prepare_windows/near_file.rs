use crate::commands::prepare_windows::{
    config::NearEdge,
    labels::{MISSING_GROUP_LABEL, validate_label_token},
};
use crate::{commands::prepare_windows::config::NearDirection, shared::io::open_text_reader};
use anyhow::{Context, Result, bail};
use fxhash::{FxHashMap, FxHashSet};
use std::{cmp::Ordering, io::BufRead, path::Path};

/// Interval from the `--near` set.
///
/// The `group` is optional; it is used when composing output groups with distance bins.
#[derive(Debug, Clone, Copy)]
pub struct NearInterval {
    pub start: u32,
    pub end: u32,
    pub group_id: Option<u32>,
    pub strand: Strand,
}

/// Per-chromosome near intervals, sorted and validated.
#[derive(Debug, Default, Clone)]
pub struct NearChrom {
    pub intervals: Vec<NearInterval>,
    pub cursor: usize,
}

/// Index holding near intervals per chromosome and a compact group-id interner.
#[derive(Debug, Default, Clone)]
pub struct NearIndex {
    pub per_chrom: FxHashMap<String, NearChrom>,
    pub group_name_to_id: FxHashMap<String, u32>,
    pub group_id_to_name: Vec<String>,
    pub warned_no_near: FxHashSet<String>,
    pub warned_no_direction: FxHashSet<String>,
}

/// Where the window sits relative to the nearest interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NearWindowSide {
    Upstream,
    Downstream,
    Overlap,
}

/// Distance, group id, and window-relative side for a selected nearest interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NearestDistance {
    pub distance: i32,
    pub group_id: Option<u32>,
    pub window_side: NearWindowSide,
}

/// Captures the left and right hits when a tie happens.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NearTie {
    pub left: Option<NearestDistance>,
    pub right: Option<NearestDistance>,
}

/// Result of the nearest-edge lookup: either one hit or a tie.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NearestResult {
    Single(NearestDistance),
    Tie(NearTie),
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum NearDuplicatesPolicy {
    /// Fail on identical (chrom,start,end) edges with a descriptive message.
    Error,
    /// Keep the first record in each run of duplicates. Drop the rest.
    KeepFirst,
    /// Drop the entire set of duplicates.
    DropAll,
    /// Merge groups across identical edges (and sometimes strands) into one record.
    ///
    /// Group names are joined with "`__`" in stable input order, with duplicates removed. Missing groups are ignored.
    Merge,
}

/// Load and validate the `--near` intervals.
///
/// Validation rules:
/// - Intervals must have end > start.
/// - Intervals must be pairwise non-overlapping per chromosome (half-open).
/// - No two intervals on the same chromosome may have identical edges.
///   When a strand column is present, identity is `(start, end, strand)`;
///   otherwise it is `(start, end)`.
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
/// - strand_col:
///     Optional column index for the strand field.
///     When omitted, all intervals default to the `+` strand.
/// - group_cols:
///     Optional column indices for group name fields.
///     When omitted, group names are left empty.
///
/// Returns
/// -------
/// - index:
///     Per-chromosome near intervals and a group interner.
pub fn load_near_index(
    path: &Path,
    separator: char,
    has_header: bool,
    strand_col: Option<usize>,
    group_cols: Option<&[usize]>,
    consider_strand_in_dups: bool,
    near_duplicates: NearDuplicatesPolicy,
) -> Result<NearIndex> {
    let mut reader = open_text_reader(path)?;
    let mut line = String::new();

    if has_header {
        line.clear();
        reader.read_line(&mut line)?; // discard header
    }

    let mut raw_by_chrom: FxHashMap<String, Vec<(u32, u32, Strand, Option<String>)>> =
        FxHashMap::default();

    for (lineno, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split(separator).collect();

        let chrom = fields
            .get(0)
            .context("Near parse: missing chrom")?
            .trim()
            .to_string();
        let start: u32 = fields
            .get(1)
            .context("Near parse: missing start")?
            .trim()
            .parse()
            .with_context(|| format!("Invalid start at near line {}", lineno + 1))?;
        let end: u32 = fields
            .get(2)
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

        let strand = if let Some(idx) = strand_col {
            let strand_value = fields
                .get(idx)
                .context("Near parse: missing strand")?
                .trim();
            match strand_value {
                "+" => Strand::Plus,
                "-" => Strand::Minus,
                "." | "" => Strand::Unknown,
                other => bail!(
                    "Near parse error at line {}: invalid strand '{}'",
                    lineno + 1,
                    other
                ),
            }
        } else {
            Strand::Plus
        };

        let group = if let Some(indices) = group_cols {
            let mut parts: Vec<&str> = Vec::with_capacity(indices.len());
            for &idx in indices {
                let name = fields.get(idx).unwrap_or(&"").trim();
                if name.is_empty() {
                    parts.push(MISSING_GROUP_LABEL);
                    continue;
                }
                if let Err(message) = validate_label_token(name, "near group label") {
                    bail!(message);
                }
                parts.push(name);
            }
            if parts.is_empty() {
                None
            } else {
                Some(parts.join("__"))
            }
        } else {
            None
        };

        raw_by_chrom
            .entry(chrom)
            .or_default()
            .push((start, end, strand, group));
    }

    // Validate and intern groups
    let mut index = NearIndex::default();

    for (chrom, mut items) in raw_by_chrom {
        // Sort by (start, end, strand, group)
        items.sort_unstable_by(|a, b| {
            a.0.cmp(&b.0)
                .then_with(|| a.1.cmp(&b.1))
                .then_with(|| a.2.cmp(&b.2))
                .then_with(|| match (&a.3, &b.3) {
                    (Some(ga), Some(gb)) => ga.cmp(gb),
                    (None, Some(_)) => Ordering::Less,
                    (Some(_), None) => Ordering::Greater,
                    (None, None) => Ordering::Equal,
                })
        });

        // Resolve identical-edge runs per policy before validating overlaps
        let items =
            resolve_identical_edges_runs(&items, near_duplicates, consider_strand_in_dups, &chrom)?;

        // Validate no overlap and no identical edges per chromosome
        let mut validated: Vec<NearInterval> = Vec::with_capacity(items.len());
        let mut last_end: u32 = 0;
        let mut have_last = false;

        for (start, end, strand, group_opt) in items.into_iter() {
            if have_last && start < last_end {
                bail!(
                    "Near validation failed on {}: intervals overlap at [{}, {}) and previous ending at {}.",
                    chrom,
                    start,
                    end,
                    last_end
                );
            }
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
                strand,
            });
        }

        index.per_chrom.insert(
            chrom,
            NearChrom {
                intervals: validated,
                cursor: 0,
            },
        );
    }

    Ok(index)
}

/// Compute distance from a window to the nearest considered `near` interval edge.
///
/// If the window overlaps the interval, the distance is 0 (we are at the site).
/// Else, the distance is the minimum of the distances from the window's two edges
/// to the selected target edges (left, right, or nearest).
///
/// Signed distances:
/// - Window upstream of the near interval (lies to the left in genomic order for `+`-strand) yields a negative distance.
/// - Window downstream of the near interval yields a positive distance.
///
/// Parameters
/// ----------
/// - `window_start`:
///     Window start.
/// - `window_end`:
///     Window end (exclusive).
/// - `near_chrom`:
///     Validated, sorted, non-overlapping intervals on the same chromosome.
/// - `which_edge`:
///     Which edge (`Left`, `Right`, `Nearest`) of the near-intervals to measure to.
/// - `directions`:
///     Whether to consider near-intervals when windows lie (`Upstream`, `Downstream`,
///     `Both`) or overlap.
/// - `signed`:
///     Whether to return signed distances.
/// - `cursor`:
///     Maintained across calls. After updating for this window, it points to the
///     last interval with `end <= window_start` if such interval exists; otherwise it can remain at 0.
///
/// Returns
/// -------
/// - result:
///     Either a single nearest hit or a tie describing left/right hits.
pub fn nearest_edge_distance(
    window_start: u32,
    window_end: u32,
    near_chrom: &mut NearChrom,
    which_edge: &NearEdge,
    directions: &NearDirection,
    signed: bool,
) -> Option<NearestResult> {
    let n = near_chrom.intervals.len();
    if n == 0 {
        return None;
    }

    // Move cursor to the last interval whose end <= window_start
    while near_chrom.cursor + 1 < n
        && near_chrom.intervals[near_chrom.cursor + 1].end <= window_start
    {
        near_chrom.cursor += 1;
    }

    // Determine left/right candidate indices based on the cursor position

    // `left_interval_idx` is index of interval left of window
    // (window is downstream for `+`-stranded interval)
    let left_interval_idx: Option<usize> =
        if near_chrom.intervals[near_chrom.cursor].end <= window_start {
            Some(near_chrom.cursor)
        } else {
            None
        };

    // `right_interval_idx` is index of interval right of window
    // (window is upstream for `+`-stranded interval)
    let right_interval_idx: Option<usize> = match left_interval_idx {
        Some(li) => {
            if li + 1 < n {
                Some(li + 1)
            } else {
                None
            }
        }
        // If None, no `left_interval_idx` existed and the cursor points to
        // the `right_interval_idx` (if such an interval exists)
        None => {
            if near_chrom.cursor < n {
                Some(near_chrom.cursor)
            } else {
                None
            }
        }
    };

    // If the right candidate overlaps the window (any overlap), we are at the site
    // By definition, the left_interval_idx cannot overlap at this point
    if let Some(ri) = right_interval_idx {
        let iv = near_chrom.intervals[ri];
        if window_start < iv.end && window_end > iv.start {
            return Some(NearestResult::Single(NearestDistance {
                distance: 0,
                group_id: iv.group_id,
                window_side: NearWindowSide::Overlap,
            }));
        }
    }

    // Helper: Evaluate the chosen edge(s) of `iv` against both window edges and
    // return the signed distance and side for the closest pairing
    #[inline]
    fn edge_distance(
        interval: &NearInterval,
        which_edge: &NearEdge,
        window_start: u32,
        window_end: u32,
    ) -> (i32, NearWindowSide) {
        // Best candidate across the considered edges of this interval
        let mut best_signed_distance: i32 = i32::MAX;
        let mut best_side: NearWindowSide = NearWindowSide::Overlap;

        // Evaluate one target edge against both window edges
        let mut consider_edge = |target_edge_bp: i32| {
            // Positive means the window lies downstream (to the right for '+' strand)
            // of the near edge. Negative means the window is upstream.
            let distance_to_window_start = window_start as i32 - target_edge_bp;
            let distance_to_window_end = window_end as i32 - target_edge_bp;

            // Choose the closer of {to-start, to-end}. Prefer start on ties
            let chosen_signed_distance =
                if distance_to_window_start.abs() <= distance_to_window_end.abs() {
                    distance_to_window_start
                } else {
                    distance_to_window_end
                };

            let side_of_target = if chosen_signed_distance < 0 {
                NearWindowSide::Upstream
            } else if chosen_signed_distance > 0 {
                NearWindowSide::Downstream
            } else {
                NearWindowSide::Overlap
            };

            if chosen_signed_distance.abs() < best_signed_distance.abs() {
                best_signed_distance = chosen_signed_distance;
                best_side = side_of_target;
            }
        };

        match which_edge {
            NearEdge::Left => consider_edge(interval.start as i32),
            NearEdge::Right => consider_edge(interval.end as i32),
            NearEdge::Nearest => {
                consider_edge(interval.start as i32);
                consider_edge(interval.end as i32);
            }
            // Use the interval edge towards upstream windows given the near interval's annotated strand orientation
            NearEdge::Upstream => match interval.strand {
                Strand::Plus => consider_edge(interval.start as i32),
                Strand::Minus => consider_edge(interval.end as i32),
                Strand::Unknown => {
                    // Fall back to genomic-nearest when strand is unknown
                    consider_edge(interval.start as i32);
                    consider_edge(interval.end as i32);
                }
            },
            // Use the interval edge towards downstream windows given the near interval's annotated strand orientation
            NearEdge::Downstream => match interval.strand {
                Strand::Plus => consider_edge(interval.end as i32),
                Strand::Minus => consider_edge(interval.start as i32),
                Strand::Unknown => {
                    // Fall back to genomic-nearest when strand is unknown
                    consider_edge(interval.start as i32);
                    consider_edge(interval.end as i32);
                }
            },
        }

        (best_signed_distance, best_side)
    }

    // Helper: Orient genomic distance/side into strand-relative space using iv.strand
    #[inline]
    fn orient_by_strand(
        genomic_side: NearWindowSide,
        genomic_dist: i32,
        iv_strand: Strand,
        signed: bool,
    ) -> (NearWindowSide, i32) {
        let (side, dist) = match iv_strand {
            Strand::Plus | Strand::Unknown => (genomic_side, genomic_dist),
            Strand::Minus => {
                let flipped_side = match genomic_side {
                    NearWindowSide::Upstream => NearWindowSide::Downstream,
                    NearWindowSide::Downstream => NearWindowSide::Upstream,
                    NearWindowSide::Overlap => NearWindowSide::Overlap,
                };
                (flipped_side, -genomic_dist)
            }
        };
        if signed {
            (side, dist)
        } else {
            (side, dist.abs())
        }
    }

    // Build at most two candidates:
    // - interval to the right of the window
    // - interval to the left of the window

    // Detect interval sitting right of the window
    let mut right_interval_hit: Option<NearestDistance> = None;
    if let Some(ri) = right_interval_idx {
        let iv = near_chrom.intervals[ri];
        let (genomic_distance, genomic_side) =
            edge_distance(&iv, which_edge, window_start, window_end);
        let (strand_relative_side, strand_relative_distance) =
            orient_by_strand(genomic_side, genomic_distance, iv.strand, signed);

        // Filter by requested direction(s). Overlap would have returned already
        let from_considered_side = match directions {
            NearDirection::Both => true,
            NearDirection::Upstream => {
                matches!(
                    strand_relative_side,
                    NearWindowSide::Upstream | NearWindowSide::Overlap
                )
            }
            NearDirection::Downstream => matches!(
                strand_relative_side,
                NearWindowSide::Downstream | NearWindowSide::Overlap
            ),
        };
        if from_considered_side {
            right_interval_hit = Some(NearestDistance {
                distance: strand_relative_distance,
                group_id: iv.group_id,
                window_side: strand_relative_side,
            });
        }
    }

    // Detect interval sitting left of the window
    let mut left_interval_hit: Option<NearestDistance> = None;
    if let Some(li) = left_interval_idx {
        let iv = near_chrom.intervals[li];
        let (genomic_distance, genomic_side) =
            edge_distance(&iv, which_edge, window_start, window_end);
        let (strand_relative_side, strand_relative_distance) =
            orient_by_strand(genomic_side, genomic_distance, iv.strand, signed);

        let from_considered_side = match directions {
            NearDirection::Both => true,
            NearDirection::Upstream => {
                matches!(
                    strand_relative_side,
                    NearWindowSide::Upstream | NearWindowSide::Overlap
                )
            }
            NearDirection::Downstream => matches!(
                strand_relative_side,
                NearWindowSide::Downstream | NearWindowSide::Overlap
            ),
        };
        if from_considered_side {
            left_interval_hit = Some(NearestDistance {
                distance: strand_relative_distance,
                group_id: iv.group_id,
                window_side: strand_relative_side,
            });
        }
    }

    // Choose winner or tie
    match (left_interval_hit, right_interval_hit) {
        (None, None) => None,
        (Some(left), None) => Some(NearestResult::Single(left)),
        (None, Some(right)) => Some(NearestResult::Single(right)),
        (Some(left), Some(right)) => {
            let left_abs_distance = left.distance.abs();
            let right_abs_distance = right.distance.abs();

            if left_abs_distance < right_abs_distance {
                Some(NearestResult::Single(left))
            } else if right_abs_distance < left_abs_distance {
                Some(NearestResult::Single(right))
            } else {
                Some(NearestResult::Tie(NearTie {
                    left: Some(left),
                    right: Some(right),
                }))
            }
        }
    }
}

#[inline]
fn resolve_identical_edges_runs(
    items_sorted: &[(u32, u32, Strand, Option<String>)],
    policy: NearDuplicatesPolicy,
    consider_strand_in_dups: bool,
    chrom: &str,
) -> Result<Vec<(u32, u32, Strand, Option<String>)>> {
    let mut out: Vec<(u32, u32, Strand, Option<String>)> = Vec::with_capacity(items_sorted.len());
    let mut i = 0usize;

    while i < items_sorted.len() {
        let t = &items_sorted[i];
        let (start, end, strand) = (t.0, t.1, t.2);
        let group_opt_ref: Option<&String> = t.3.as_ref();

        // Advance j over the run of identical keys
        let mut j = i + 1;
        while j < items_sorted.len() {
            let t2 = &items_sorted[j];
            let same_edge = t2.0 == start && t2.1 == end;
            let same_key = if consider_strand_in_dups {
                same_edge && t2.2 == strand
            } else {
                same_edge
            };
            if same_key {
                j += 1;
            } else {
                break;
            }
        }

        let run_len = j - i;
        if run_len == 1 {
            // No duplicates -> keep as-is
            out.push((start, end, strand, group_opt_ref.cloned()));
        } else {
            match policy {
                NearDuplicatesPolicy::Error => {
                    let mut groups: Vec<String> = items_sorted[i..j]
                        .iter()
                        .filter_map(|t| t.3.as_ref().cloned())
                        .collect();
                    groups.sort();
                    groups.dedup();
                    let groups_display = if groups.is_empty() {
                        ".".to_string()
                    } else {
                        groups.join(", ")
                    };
                    bail!(
                        "Near validation failed on {}: duplicate edges at [{}, {}) (strand {:?}). Found groups: {}. Use --near-duplicates to resolve.",
                        chrom,
                        start,
                        end,
                        strand,
                        groups_display
                    );
                }
                NearDuplicatesPolicy::KeepFirst => {
                    out.push((start, end, strand, group_opt_ref.cloned()));
                }
                NearDuplicatesPolicy::DropAll => {
                    // Keep none from this run
                }
                NearDuplicatesPolicy::Merge => {
                    use fxhash::FxHashSet;
                    let mut seen: FxHashSet<&str> = FxHashSet::default();
                    let mut merged: Vec<&str> = Vec::new();
                    for t in &items_sorted[i..j] {
                        if let Some(name) = t.3.as_deref() {
                            if seen.insert(name) {
                                merged.push(name);
                            }
                        }
                    }
                    let merged_opt = if merged.is_empty() {
                        None
                    } else {
                        Some(merged.join("__"))
                    };
                    out.push((start, end, strand, merged_opt));
                }
            }
        }

        i = j;
    }

    Ok(out)
}
