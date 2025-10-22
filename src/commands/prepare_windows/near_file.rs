use crate::commands::prepare_windows::config::NearEdge;
use crate::{commands::prepare_windows::config::NearDirection, shared::io::open_text_reader};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
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
#[derive(Debug, Default)]
pub struct NearChrom {
    pub intervals: Vec<NearInterval>,
    pub cursor: usize,
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

// TODO: Rename ("hit" is a bad name)
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
/// - strand_col_present:
///     Whether the strand column is present (fourth column).
/// - group_col_present:
///     Whether a `group` column is present (fourth or fifth column).
///
/// Returns
/// -------
/// - index:
///     Per-chromosome near intervals and a group interner.
pub fn load_near_index(
    path: &Path,
    separator: char,
    has_header: bool,
    strand_col_present: bool,
    group_col_present: bool,
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

        let strand = if strand_col_present {
            match it.next().map(|s| s.trim()) {
                Some("+") => Strand::Plus,
                Some("-") => Strand::Minus,
                Some(".") | None => Strand::Unknown,
                Some(other) => bail!(
                    "Near parse error at line {}: invalid strand '{}'",
                    lineno + 1,
                    other
                ),
            }
        } else {
            Strand::Plus // Default when not provided
        };

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
/// - Target interval being *upstream* of the window gives a negative distance. *Downstream* gives a positive distance.
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
///     Which direction(s) (`Upstream`, `Downstream`, `Both`) to consider near-intervals at.
///     E.g., `Upstream` only considers near-intervals that lie before (or overlap) the window.
/// - `signed`:
///     Whether to return signed distances.
/// - `cursor`:
///     Maintained across calls. After updating for this window, it points to the
///     last interval with `end <= window_start` if such interval exists; otherwise it can remain at 0.
///
/// Returns
/// -------
/// - result:
///     Either a single nearest hit or a tie describing upstream/downstream hits.
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

    // Move cursor to the last interval whose end <= window_start.
    while near_chrom.cursor + 1 < n
        && near_chrom.intervals[near_chrom.cursor + 1].end <= window_start
    {
        near_chrom.cursor += 1;
    }

    // Determine upstream/downstream candidate indices based on the cursor position.
    let upstream_idx: Option<usize> = if near_chrom.intervals[near_chrom.cursor].end <= window_start
    {
        Some(near_chrom.cursor)
    } else {
        None
    };
    let downstream_idx: Option<usize> = match upstream_idx {
        Some(ui) => {
            if ui + 1 < n {
                Some(ui + 1)
            } else {
                None
            }
        }
        None => {
            if near_chrom.cursor < n {
                Some(near_chrom.cursor)
            } else {
                None
            }
        }
    };

    // If the downstream candidate overlaps the window (any overlap), we are at the site
    if let Some(di) = downstream_idx {
        let iv = near_chrom.intervals[di];
        if window_start < iv.end && window_end > iv.start {
            return Some(NearestResult::Single(NearHit {
                distance: 0,
                group_id: iv.group_id,
                side: NearSide::Overlap,
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
    ) -> (i32, NearSide) {
        // Best candidate across the considered edges of this interval
        let mut best_signed_distance: i32 = i32::MAX;
        let mut best_side: NearSide = NearSide::Overlap;

        // Evaluate one target edge against both window edges
        let mut consider_edge = |target_edge_bp: i32| {
            let distance_to_window_start = target_edge_bp - window_start as i32;
            let distance_to_window_end = target_edge_bp - window_end as i32;

            // Choose the closer of {to-start, to-end}; prefer start on ties
            let chosen_signed_distance =
                if distance_to_window_start.abs() <= distance_to_window_end.abs() {
                    distance_to_window_start
                } else {
                    distance_to_window_end
                };

            let side_of_target = if chosen_signed_distance < 0 {
                NearSide::Upstream
            } else if chosen_signed_distance > 0 {
                NearSide::Downstream
            } else {
                NearSide::Overlap
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
            // Use the edge that is upstream of the near interval given its annotated strand orientation.
            NearEdge::Upstream => match interval.strand {
                Strand::Plus => consider_edge(interval.start as i32),
                Strand::Minus => consider_edge(interval.end as i32),
                Strand::Unknown => {
                    // Fall back to genomic-nearest when strand is unknown
                    consider_edge(interval.start as i32);
                    consider_edge(interval.end as i32);
                }
            },
            // Use the edge that is downstream of the near interval given its annotated strand orientation.
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
        genomic_side: NearSide,
        genomic_dist: i32,
        iv_strand: Strand,
        signed: bool,
    ) -> (NearSide, i32) {
        let (side, dist) = match iv_strand {
            Strand::Plus | Strand::Unknown => (genomic_side, genomic_dist),
            Strand::Minus => {
                let flipped_side = match genomic_side {
                    NearSide::Upstream => NearSide::Downstream,
                    NearSide::Downstream => NearSide::Upstream,
                    NearSide::Overlap => NearSide::Overlap,
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

    // Build at most two candidates: upstream and downstream.
    let mut upstream_hit: Option<NearHit> = None;
    if let Some(ui) = upstream_idx {
        let iv = near_chrom.intervals[ui];
        let (genomic_distance, genomic_side) =
            edge_distance(&iv, which_edge, window_start, window_end);
        let (strand_relative_side, strand_relative_distance) =
            orient_by_strand(genomic_side, genomic_distance, iv.strand, signed);

        // Filter by requested direction(s). Overlap would have returned already.
        let from_considered_side = match directions {
            NearDirection::Both => true,
            NearDirection::Upstream => {
                matches!(strand_relative_side, NearSide::Upstream | NearSide::Overlap)
            }
            NearDirection::Downstream => matches!(
                strand_relative_side,
                NearSide::Downstream | NearSide::Overlap
            ),
        };
        if from_considered_side {
            upstream_hit = Some(NearHit {
                distance: strand_relative_distance,
                group_id: iv.group_id,
                side: strand_relative_side,
            });
        }
    }

    let mut downstream_hit: Option<NearHit> = None;
    if let Some(di) = downstream_idx {
        let iv = near_chrom.intervals[di];
        let (genomic_distance, genomic_side) =
            edge_distance(&iv, which_edge, window_start, window_end);
        let (strand_relative_side, strand_relative_distance) =
            orient_by_strand(genomic_side, genomic_distance, iv.strand, signed);

        let from_considered_side = match directions {
            NearDirection::Both => true,
            NearDirection::Upstream => {
                matches!(strand_relative_side, NearSide::Upstream | NearSide::Overlap)
            }
            NearDirection::Downstream => matches!(
                strand_relative_side,
                NearSide::Downstream | NearSide::Overlap
            ),
        };
        if from_considered_side {
            downstream_hit = Some(NearHit {
                distance: strand_relative_distance,
                group_id: iv.group_id,
                side: strand_relative_side,
            });
        }
    }

    // Choose winner or tie
    match (upstream_hit, downstream_hit) {
        (None, None) => None,
        (Some(upstream), None) => Some(NearestResult::Single(upstream)),
        (None, Some(downstream)) => Some(NearestResult::Single(downstream)),
        (Some(upstream), Some(downstream)) => {
            let upstream_abs_distance = upstream.distance.abs();
            let downstream_abs_distance = downstream.distance.abs();

            if upstream_abs_distance < downstream_abs_distance {
                Some(NearestResult::Single(upstream))
            } else if downstream_abs_distance < upstream_abs_distance {
                Some(NearestResult::Single(downstream))
            } else {
                Some(NearestResult::Tie(NearTie {
                    upstream: Some(upstream),
                    downstream: Some(downstream),
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
