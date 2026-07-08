use crate::{
    commands::cli_common::{
        DistributionWindowSpec, WindowAssigner, min_overlap_fraction_for_window_assignment,
    },
    shared::{
        base::ZEROISH_F64_TOLERANCE,
        interval::Interval,
        kmers::{
            kmer_codec::{Kmer, KmerCodes, KmerOrientation},
            motifs_file::{EncodedMotifKey, SelectedMotifLookup},
        },
        overlaps::{FixedWidthOverlapCursor, FixedWidthWindowSource, TileBedWindowView},
        windowing::DistributionWindowContext,
    },
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use std::ops::Range;

/// Sparse full-space reference k-mer counts for one output row.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct KmerCounts {
    pub(crate) counts: FxHashMap<Kmer, f64>,
}

impl KmerCounts {
    /// Increment one encoded k-mer by a count weight.
    #[inline]
    pub(crate) fn incr_weighted(&mut self, kmer: Kmer, weight: f64) {
        *self.counts.entry(kmer).or_insert(0.0) += weight;
    }

    /// Return whether a weight should create a sparse count entry.
    #[inline]
    pub(crate) fn should_store_weight(weight: f64) -> Result<bool> {
        ensure!(
            weight.is_finite(),
            "sparse reference k-mer weight {weight} is not finite"
        );
        ensure!(
            weight >= -ZEROISH_F64_TOLERANCE,
            "sparse reference k-mer weight {weight} is negative, this is not currently supported"
        );
        Ok(weight > ZEROISH_F64_TOLERANCE)
    }
}

/// Sparse full-space reference k-mer counts keyed by global output row.
pub(crate) type KmerCountsByWindow = FxHashMap<u64, KmerCounts>;

/// Sparse motifs-file target counts for one output row.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct SelectedKmerCounts {
    pub(crate) counts: FxHashMap<u32, f64>,
}

impl SelectedKmerCounts {
    /// Increment one selected motifs-file target by a count weight.
    #[inline]
    pub(crate) fn incr_weighted(&mut self, target_idx: u32, weight: f64) {
        *self.counts.entry(target_idx).or_insert(0.0) += weight;
    }
}

impl IntoIterator for SelectedKmerCounts {
    type Item = (u32, f64);
    type IntoIter = std::collections::hash_map::IntoIter<u32, f64>;

    fn into_iter(self) -> Self::IntoIter {
        self.counts.into_iter()
    }
}

/// Sparse motifs-file target counts keyed by global output row.
pub(crate) type SelectedKmerCountsByWindow = FxHashMap<u64, SelectedKmerCounts>;

/// Count k-mers owned by one tile into every overlapping output window.
///
/// `owned_starts` is the range of k-mer start offsets in `enc.codes` that this tile owns, relative
/// to the loaded reference sequence whose chromosome start is `sequence_start`. Each k-mer interval
/// must have a non-sentinel code and fit completely within `chrom_len` before window lookup.
/// `window_source` must already be prepared for this tile. Fixed-size windows use
/// `window_context` to add chromosome row offsets, while BED-like windows must carry their output
/// row ids in the overlap records produced by the prepared source.
pub(crate) fn count_kmers_by_window(
    counts_by_window: &mut KmerCountsByWindow,
    selected_counts_by_window: &mut SelectedKmerCountsByWindow,
    enc: &Enc<'_>,
    window_context: &DistributionWindowContext<'_>,
    window_source: FixedWidthWindowSource,
    owned_starts: Range<u64>,
    sequence_start: u64,
    chrom_len: u64,
    assign_by: WindowAssigner,
    selected_motifs: Option<&SelectedMotifLookup>,
) -> Result<()> {
    let k = enc.k as u64;
    let min_overlap_fraction = min_overlap_fraction_for_window_assignment(assign_by, k);
    let query_width = query_width_for_assignment(assign_by, k);
    let mut overlap_cursor =
        FixedWidthOverlapCursor::new(chrom_len, window_source, query_width, min_overlap_fraction)?;

    for kmer_start in owned_starts {
        let code = enc.codes.get(kmer_start as usize);
        if code == enc.none || code == enc.n {
            continue;
        }

        let kmer_start_abs = sequence_start
            .checked_add(kmer_start)
            .context("reference k-mer start coordinate overflowed")?;
        let kmer_end_abs = kmer_start_abs
            .checked_add(k)
            .context("reference k-mer end coordinate overflowed")?;
        if kmer_end_abs > chrom_len {
            continue;
        }

        let selected_target_idx = if let Some(selected_motifs) = selected_motifs {
            let key = EncodedMotifKey {
                inside_code: code,
                outside_code: 0,
                reverse_on_decode: false,
            };
            match selected_motifs.target_for(key) {
                Some(target_idx) => Some(target_idx),
                None => continue,
            }
        } else {
            None
        };

        let query_start = window_query_start(assign_by, kmer_start_abs, k)?;
        let Some(overlapping_windows) = overlap_cursor.find_overlaps(query_start)? else {
            continue;
        };

        for overlapped_window in overlapping_windows.windows {
            let row_idx = match window_context.spec {
                DistributionWindowSpec::Bed(_) | DistributionWindowSpec::GroupedBed(_) => {
                    overlapped_window
                        .output_idx
                        .context("BED reference k-mer overlap did not carry an output row id")?
                }
                DistributionWindowSpec::Global | DistributionWindowSpec::Size(_) => {
                    window_context.original_idx(overlapped_window.idx)
                }
            };
            let weight = count_weight(assign_by, overlapped_window.overlap_fraction);
            if let Some(target_idx) = selected_target_idx {
                selected_counts_by_window
                    .entry(row_idx)
                    .or_default()
                    .incr_weighted(target_idx, weight);
            } else {
                counts_by_window
                    .entry(row_idx)
                    .or_default()
                    .incr_weighted(forward_kmer(enc.k, code), weight);
            }
        }
    }
    Ok(())
}

/// Prepare the fixed-width window source for one tile before k-mer counting starts.
///
/// Global and fixed-size modes can be represented directly. BED-like modes use the tile-local BED
/// window view precomputed by the command runner so the overlap cursor only sees the windows that
/// can matter for this tile. A missing BED view means the chromosome has no BED windows.
///
/// `owned_starts` uses offsets into the loaded tile sequence, not absolute chromosome
/// coordinates. The helper clips that range to k-mer starts whose full `[start, start + k)` span
/// fits inside `chrom_len`, then converts the first and last surviving starts into an absolute
/// assignment envelope for the BED source. That envelope is what lets the overlap cursor separate
/// rows that always cover this tile's queries from rows that still need per-query checks.
pub(crate) fn prepare_ref_kmer_window_source(
    window_context: &DistributionWindowContext<'_>,
    tile_bed_window_view: Option<TileBedWindowView<'_>>,
    owned_starts: Range<u64>,
    sequence_start: u64,
    chrom_len: u64,
    kmer_size: u64,
    assign_by: WindowAssigner,
) -> Result<Option<FixedWidthWindowSource>> {
    match window_context.spec {
        DistributionWindowSpec::Global => Ok(Some(FixedWidthWindowSource::Global)),
        DistributionWindowSpec::Size(window_bp) => {
            Ok(Some(FixedWidthWindowSource::FixedSize(*window_bp)))
        }
        DistributionWindowSpec::Bed(_) | DistributionWindowSpec::GroupedBed(_) => {
            let Some(tile_bed_window_view) = tile_bed_window_view else {
                return Ok(None);
            };
            if tile_bed_window_view
                .chromosome_windows
                .all_windows
                .is_empty()
            {
                return Ok(None);
            }
            let query_width = query_width_for_assignment(assign_by, kmer_size);
            let Some(valid_owned_starts) = valid_owned_starts_for_full_kmers(
                owned_starts,
                sequence_start,
                chrom_len,
                kmer_size,
            )?
            else {
                return Ok(None);
            };
            let tile_assignment_envelope = tile_assignment_envelope(
                assign_by,
                kmer_size,
                query_width,
                sequence_start,
                valid_owned_starts,
            )?;

            Ok(Some(FixedWidthWindowSource::bed_from_tile_view(
                chrom_len,
                tile_bed_window_view,
                tile_assignment_envelope,
            )?))
        }
    }
}

/// Return the fixed query width used by the overlap cursor for one assignment mode.
///
/// Midpoint assignment projects each k-mer to a single center base. Other modes compare the full
/// k-mer span against candidate windows.
fn query_width_for_assignment(assign_by: WindowAssigner, kmer_size: u64) -> u64 {
    match assign_by {
        WindowAssigner::Midpoint => 1,
        WindowAssigner::CountOverlap
        | WindowAssigner::Any
        | WindowAssigner::All
        | WindowAssigner::Proportion(_) => kmer_size,
    }
}

/// Return the tile-owned starts that can form complete k-mers on this chromosome.
///
/// Input and output ranges are relative to the loaded reference sequence for the tile. The
/// chromosome coordinate of a start is `sequence_start + start`. Starts whose complete k-mer would
/// extend beyond `chrom_len` are removed here, matching the counting loop's later safety check.
///
/// A start `s` is usable when the complete k-mer fits inside the chromosome:
///
/// `sequence_start + s + kmer_size <= chrom_len`
///
/// Equivalently, when `kmer_size <= chrom_len`, the last valid absolute start is:
///
/// `last_valid_abs = chrom_len - kmer_size`
///
/// If the first owned start is after `last_valid_abs`, no owned start can form a full k-mer.
/// Otherwise, clip the exclusive end of `owned_starts` to one past the last valid relative start:
///
/// `end = min(owned_starts.end, last_valid_abs - sequence_start + 1)`
///
/// and return `owned_starts.start..end`. The implementation computes the same clipping through
/// absolute coordinates and uses checked arithmetic so overflow or inconsistent coordinates fail
/// instead of wrapping.
fn valid_owned_starts_for_full_kmers(
    owned_starts: Range<u64>,
    sequence_start: u64,
    chrom_len: u64,
    kmer_size: u64,
) -> Result<Option<Range<u64>>> {
    if owned_starts.is_empty() || kmer_size > chrom_len {
        return Ok(None);
    }

    let max_valid_start_abs = chrom_len - kmer_size;
    let first_start_abs = sequence_start
        .checked_add(owned_starts.start)
        .context("first owned reference k-mer start coordinate overflowed")?;
    if first_start_abs > max_valid_start_abs {
        return Ok(None);
    }

    let last_owned_start = owned_starts
        .end
        .checked_sub(1)
        .context("owned reference k-mer start range was empty")?;
    let last_start_abs = sequence_start
        .checked_add(last_owned_start)
        .context("last owned reference k-mer start coordinate overflowed")?
        .min(max_valid_start_abs);
    let last_valid_start = last_start_abs
        .checked_sub(sequence_start)
        .context("last valid reference k-mer start fell before the loaded sequence")?;

    Ok(Some(owned_starts.start..last_valid_start + 1))
}

/// Build the absolute assignment envelope for every valid k-mer start in a tile.
///
/// The returned interval is in chromosome coordinates and contains every query interval that
/// `count_kmers_by_window` can submit to the overlap cursor for this tile. For midpoint mode, the
/// envelope spans the first through last midpoint base. For all other assignment modes, it spans
/// the first through last complete k-mer interval.
fn tile_assignment_envelope(
    assign_by: WindowAssigner,
    kmer_size: u64,
    query_width: u64,
    sequence_start: u64,
    valid_owned_starts: Range<u64>,
) -> Result<Interval<u64>> {
    let first_start_abs = sequence_start
        .checked_add(valid_owned_starts.start)
        .context("first reference k-mer assignment coordinate overflowed")?;
    let last_start_abs = sequence_start
        .checked_add(valid_owned_starts.end - 1)
        .context("last reference k-mer assignment coordinate overflowed")?;
    let first_query_start = window_query_start(assign_by, first_start_abs, kmer_size)?;
    let last_query_start = window_query_start(assign_by, last_start_abs, kmer_size)?;
    let envelope_end = last_query_start
        .checked_add(query_width)
        .context("reference k-mer assignment envelope end overflowed")?;

    Ok(Interval::new(first_query_start, envelope_end)?)
}

fn window_query_start(assign_by: WindowAssigner, kmer_start: u64, kmer_size: u64) -> Result<u64> {
    match assign_by {
        WindowAssigner::Midpoint => {
            ensure!(
                kmer_size % 2 == 1,
                "`--assign-by midpoint` requires an odd `--kmer-size`"
            );
            kmer_start
                .checked_add(kmer_size / 2)
                .context("k-mer midpoint coordinate overflowed")
        }
        WindowAssigner::CountOverlap
        | WindowAssigner::Any
        | WindowAssigner::All
        | WindowAssigner::Proportion(_) => Ok(kmer_start),
    }
}

fn count_weight(assign_by: WindowAssigner, overlap_fraction: f64) -> f64 {
    match assign_by {
        WindowAssigner::CountOverlap => overlap_fraction,
        WindowAssigner::Any
        | WindowAssigner::All
        | WindowAssigner::Midpoint
        | WindowAssigner::Proportion(_) => 1.0,
    }
}

#[inline]
fn forward_kmer(k: u8, code: u64) -> Kmer {
    Kmer {
        k,
        code,
        orientation: KmerOrientation::Forward,
    }
}

/// Container for storing k, codes, and sentinels
pub(crate) struct Enc<'a> {
    pub(crate) k: u8,
    pub(crate) codes: &'a KmerCodes,
    pub(crate) none: u64,
    pub(crate) n: u64,
}

#[cfg(test)]
mod tests {
    include!("counting_tests.rs");
}
