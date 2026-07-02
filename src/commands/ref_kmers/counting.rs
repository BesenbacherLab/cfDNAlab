use crate::{
    commands::cli_common::{WindowAssigner, min_overlap_fraction_for_window_assignment},
    shared::{
        base::ZEROISH_F64_TOLERANCE,
        kmers::{
            kmer_codec::{Kmer, KmerCodes, KmerOrientation},
            motifs_file::{EncodedMotifKey, SelectedMotifLookup},
        },
        overlaps::FixedWidthOverlapCursor,
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
/// must have a non-sentinel code and fit completely within `chrom_len` before window lookup. Output
/// row mapping goes through `window_context`, so fixed windows use chromosome offsets and BED-like
/// windows use their stored row ids.
pub(crate) fn count_kmers_by_window(
    counts_by_window: &mut KmerCountsByWindow,
    selected_counts_by_window: &mut SelectedKmerCountsByWindow,
    enc: &Enc<'_>,
    window_context: &DistributionWindowContext<'_>,
    window_pointer: &mut usize,
    owned_starts: Range<u64>,
    sequence_start: u64,
    chrom_len: u64,
    assign_by: WindowAssigner,
    selected_motifs: Option<&SelectedMotifLookup>,
) -> Result<()> {
    let k = enc.k as u64;
    let min_overlap_fraction = min_overlap_fraction_for_window_assignment(assign_by, k);
    let query_width = match assign_by {
        WindowAssigner::Midpoint => 1,
        WindowAssigner::CountOverlap
        | WindowAssigner::Any
        | WindowAssigner::All
        | WindowAssigner::Proportion(_) => k,
    };
    if window_context.requires_windows()
        && window_context
            .windows_slice()
            .map_or(true, |windows| windows.is_empty())
    {
        return Ok(());
    }
    let mut overlap_cursor = FixedWidthOverlapCursor::new(
        chrom_len,
        window_context.windows_slice(),
        window_context.by_size(),
        query_width,
        min_overlap_fraction,
        *window_pointer,
    )?;

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
            let row_idx = window_context.original_idx(overlapped_window.idx);
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
