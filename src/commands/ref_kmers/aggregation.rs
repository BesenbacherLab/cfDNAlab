use crate::{
    commands::ref_kmers::{
        counting::{KmerCounts, KmerCountsByWindow, SelectedKmerCountsByWindow},
        tiling::{
            TileResult, deserialize_selected_tile_counts, deserialize_tile_counts,
            merge_selected_tile_count_records, merge_tile_count_records,
        },
        zarr::{
            RefKmerFrequencyBins, normalize_count_bins_to_frequencies, postprocess_ref_kmer_counts,
        },
    },
    shared::kmers::{
        kmer_codec::KmerSpec,
        motifs_file::{SelectedMotifColumnKind, SelectedMotifLookup},
        process_counts::postprocess_selected_motif_counts,
    },
};
use anyhow::{Context, Result};

pub(crate) struct CollectedRefKmerFrequencies {
    pub(crate) frequency_bins: RefKmerFrequencyBins,
    pub(crate) motif_order: Vec<String>,
    pub(crate) column_kind: SelectedMotifColumnKind,
}

/// Collect tile count files and convert them into the final reference k-mer frequency matrix.
///
/// This is the join point between tile-level counting and Zarr writing. Full-space output keeps
/// encoded k-mer keys until `postprocess_ref_kmer_counts` decodes and optionally canonicalizes
/// them. Motifs-file output has already reduced counts to parser-assigned target indices, so it
/// uses the shared selected-motif postprocessor and then applies the same row normalization as the
/// full-space path.
///
/// Parameters
/// ----------
/// - `tile_results`:
///   Per-tile count files produced by the counting step
/// - `selected_motifs`:
///   Optional motifs-file lookup. When present, counts are target-indexed rather than k-mer-keyed
/// - `kmer_decode_spec`:
///   Codec for full-space output. This is required only when `selected_motifs` is absent
/// - `total_windows`:
///   Number of rows in the final output matrix
/// - `canonical`:
///   Whether full-space motif labels should be reverse-complement collapsed
/// - `all_motifs`:
///   Whether to retain the full motif axis or all motifs-file targets
///
/// Returns
/// -------
/// - `Result<CollectedRefKmerFrequencies>`:
///   Frequency rows, motif labels, and the kind of motif axis to write
pub(crate) fn collect_ref_kmer_frequencies(
    tile_results: &[TileResult],
    selected_motifs: Option<&SelectedMotifLookup>,
    kmer_decode_spec: Option<&KmerSpec>,
    total_windows: usize,
    canonical: bool,
    all_motifs: bool,
) -> Result<CollectedRefKmerFrequencies> {
    match selected_motifs {
        None => {
            let kmer_decode_spec = kmer_decode_spec
                .as_ref()
                .context("missing k-mer decode spec for full reference k-mer output")?;
            let mut reduced_counts = KmerCountsByWindow::default();
            for tile_result in tile_results {
                let count_records = deserialize_tile_counts(&tile_result.counts_path)?;
                merge_tile_count_records(&mut reduced_counts, count_records)?;
            }
            let (frequency_bins, motif_order) = postprocess_ref_kmer_counts(
                reduced_counts,
                total_windows,
                kmer_decode_spec,
                canonical,
                all_motifs,
            )?;
            Ok(CollectedRefKmerFrequencies {
                frequency_bins,
                motif_order,
                column_kind: SelectedMotifColumnKind::Motif,
            })
        }
        Some(lookup) => {
            let mut reduced_selected_counts = SelectedKmerCountsByWindow::default();
            for tile_result in tile_results {
                let count_records = deserialize_selected_tile_counts(&tile_result.counts_path)?;
                merge_selected_tile_count_records(&mut reduced_selected_counts, count_records)?;
            }
            let (count_bins, motif_order) = postprocess_selected_motif_counts(
                reduced_selected_counts,
                total_windows,
                &lookup.labels,
                all_motifs,
                KmerCounts::should_store_weight,
                crate::commands::ref_kmers::zarr::ensure_dense_ref_kmer_output_size,
            )?;
            let frequency_bins = normalize_count_bins_to_frequencies(count_bins)?;
            Ok(CollectedRefKmerFrequencies {
                frequency_bins,
                motif_order,
                column_kind: lookup.column_kind,
            })
        }
    }
}
