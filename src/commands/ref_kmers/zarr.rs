//! Reference k-mer Zarr writer.
//!
//! The public values are frequencies, not counts:
//!
//! - dense output writes `frequencies[row, motif]`
//! - sparse output writes `sparse/{row,motif,frequency,shape,sparse_dimension}`
//! - `row_scaling_factor[row]` stores the row total needed to reconstruct counts as
//!   `frequency * row_scaling_factor`
//!
//! Row and motif metadata are stored in the package so downstream code can interpret the matrix
//! without joining sidecar TSV files.

use crate::{
    commands::{
        cli_common::WindowAssigner,
        ref_kmers::counting::{KmerCounts, KmerCountsByWindow},
    },
    shared::{
        base::make_canonical,
        bed::GroupedWindows,
        blacklist::compute_blacklist_overlap,
        interval::Interval,
        io::dot_join,
        kmers::{
            kmer_codec::{Kmer, KmerOrientation, KmerSpec},
            motifs_file::SelectedMotifColumnKind,
            process_counts::{acgt_radix5_code_from_radix4, all_motifs as all_ref_kmer_motifs},
        },
        reference::ContigFootprintEntry,
        windowing::WindowBinInfo,
        zarr::{
            ZARR_ASCII_FILL_VALUE, ZARR_FLOAT64_FILL_VALUE, ZARR_INT32_FILL_VALUE,
            ZARR_INT64_FILL_VALUE, checked_i32, checked_i64, checked_index_axis, create_zarr_array,
            create_zarr_store, validate_zarr_label, write_single_chunk_zarr_array,
            write_zarr_group_metadata, write_zarr_root_metadata,
        },
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use ndarray::{Array2, ArrayView2};
use serde_json::json;
use std::{
    collections::BTreeSet,
    env,
    path::{Path, PathBuf},
    sync::Arc,
};
use zarrs::{
    array::{Array, ArraySubset, data_type},
    filesystem::FilesystemStore,
};

const CFDNALAB_REF_KMER_SCHEMA_VERSION: u32 = 1;
const TARGET_DENSE_FREQUENCY_CHUNK_CELLS: usize = 2_000_000;
/// Target maximum number of COO entries written per sparse Zarr chunk.
///
/// The sparse writer stores `row`, `motif`, and `frequency` as aligned one-dimensional arrays.
/// Keeping this limit below the dense chunk target bounds temporary memory while still writing
/// large contiguous slices to the store.
const TARGET_SPARSE_COO_CHUNK_ENTRIES: usize = 1_000_000;
const DEFAULT_MAX_DENSE_REF_KMER_OUTPUT_BYTES: u64 = 5 * 1024 * 1024 * 1024;
const MAX_DENSE_REF_KMER_OUTPUT_BYTES_ENV: &str = "CFDNALAB_REF_KMERS_MAX_DENSE_OUTPUT_BYTES";

/// Sparse count row keyed by final output column index.
pub(crate) type RefKmerCountBin = FxHashMap<u32, f64>;

/// Sparse frequency row keyed by final output column index.
pub(crate) type RefKmerFrequencyBin = FxHashMap<u32, f64>;

/// Normalized reference k-mer matrix values and their row count totals.
#[derive(Debug)]
pub(crate) struct RefKmerFrequencyBins {
    /// Per-row motif frequencies keyed by final motif-column index
    pub(crate) frequency_bins: Vec<RefKmerFrequencyBin>,
    /// Per-row total count used to reconstruct counts
    pub(crate) row_scaling_factors: Vec<f64>,
}

/// Row metadata for the public `row` axis.
pub(crate) enum RefKmerRowMetadata<'a> {
    /// One row covering all selected chromosomes.
    Global,
    /// Genomic window rows in count-row order.
    Windows {
        /// Per-row chromosome, coordinate, and blacklist metadata
        bin_info: &'a [WindowBinInfo],
        /// Whether these rows came from size windows or BED windows
        row_mode: RefKmerWindowRowMode,
    },
    /// Grouped BED rows in count-row order.
    Groups(Vec<RefKmerGroupSummary<'a>>),
}

/// Source of genomic window count rows.
pub(crate) enum RefKmerWindowRowMode {
    /// Generated fixed-size windows
    Size,
    /// User-provided BED windows
    Bed,
}

/// One grouped-BED output row.
#[derive(Debug)]
pub(crate) struct RefKmerGroupSummary<'a> {
    /// Zero-based group index, expected to match the count row index
    pub(crate) group_idx: u64,
    /// Public group label from the grouped BED input
    pub(crate) group_name: &'a str,
    /// Number of grouped BED windows contributing to the group
    pub(crate) eligible_windows: usize,
    /// Fraction of blacklisted bases across this group's windows.
    ///
    /// Groups with zero total bases use `0.0`.
    pub(crate) blacklisted_fraction: f64,
}

/// Inputs needed to write a reference k-mer frequency Zarr package.
pub(crate) struct RefKmerZarrPackage<'a> {
    /// Per-row motif frequencies keyed by final motif-column index.
    ///
    /// These values are expected to have been normalized upstream, usually by
    /// `normalize_count_bins_to_frequencies`, so each stored value is a finite non-negative row
    /// frequency and the corresponding row total is stored in `row_scaling_factors`.
    pub(crate) frequency_bins: &'a [RefKmerFrequencyBin],
    /// Per-row count totals used to reconstruct counts from frequencies.
    pub(crate) row_scaling_factors: &'a [f64],
    pub(crate) motif_labels: &'a [String],
    pub(crate) column_kind: SelectedMotifColumnKind,
    pub(crate) row_metadata: RefKmerRowMetadata<'a>,
    /// Whether to write the rectangular `frequencies[row, motif]` array
    ///
    /// This can be true even when `all_motifs` is false, if the motif axis is already complete.
    pub(crate) write_dense_output: bool,
    pub(crate) kmer_size: u8,
    pub(crate) canonical: bool,
    /// Whether the user requested `--all-motifs`, stored as output metadata
    pub(crate) all_motifs: bool,
    pub(crate) assign_by: WindowAssigner,
    pub(crate) reference_contig_footprint: &'a [ContigFootprintEntry],
}

/// Decode full reference k-mer counts into the public motif matrix.
///
/// Tile counting stores compact `Kmer` keys and sparse row counts. This function is the boundary
/// where those internal keys become output motif labels, then output column indices, then
/// row-normalized frequencies. It preserves the final output row IDs assigned before tiling.
///
/// The motif axis has two modes. With `all_motifs`, the axis contains every A/C/G/T k-mer for the
/// requested k and missing motifs are represented as zeroes in dense output. Without
/// `all_motifs`, the axis contains only motifs that were observed after filtering and optional
/// canonicalization, sorted deterministically for stable output.
///
/// Decoded motifs containing `N` are skipped before normalization because the public reference
/// k-mer axis is defined over concrete A/C/G/T labels. When canonical output is requested, labels
/// are canonicalized before aggregation so reverse-complement partners contribute to the same
/// output column.
///
/// Parameters
/// ----------
/// - `counts_by_window`:
///   Sparse tile-reduced counts keyed by final output row and encoded reference k-mer
/// - `total_windows`:
///   Number of rows in the final output matrix, including rows with no retained motifs
/// - `kmer_spec`:
///   Codec used to decode the compact k-mer keys into motif strings
/// - `canonical`:
///   Whether to collapse each motif with its reverse complement before output
/// - `all_motifs`:
///   Whether to retain all A/C/G/T k-mers instead of only observed motifs
///
/// Returns
/// -------
/// - `Result<(RefKmerFrequencyBins, Vec<String>)>`:
///   Row-normalized frequencies keyed by output column, plus the motif labels for those columns
pub(crate) fn postprocess_ref_kmer_counts(
    counts_by_window: KmerCountsByWindow,
    total_windows: usize,
    kmer_spec: &KmerSpec,
    canonical: bool,
    all_motifs: bool,
) -> Result<(RefKmerFrequencyBins, Vec<String>)> {
    // Keep string-keyed bins while decoding, filtering, and canonicalizing labels
    let mut label_bins = vec![FxHashMap::default(); total_windows];

    for (original_idx, counts) in counts_by_window {
        // Validate row IDs before indexing because tile results are deserialized inputs here
        let row_idx: usize = original_idx
            .try_into()
            .context("reference k-mer window index does not fit in usize")?;
        ensure!(
            row_idx < label_bins.len(),
            "reference k-mer window index {} is out of bounds for {} output windows",
            row_idx,
            label_bins.len()
        );

        for (kmer, value) in counts.counts {
            // Drop zeroish weights, reject NaN or invalid weights, and keep row totals meaningful
            if !KmerCounts::should_store_weight(value)? {
                continue;
            }
            ensure!(
                usize::from(kmer.k) == kmer_spec.k,
                "reference k-mer key has k={} but output k-mer size is {}",
                kmer.k,
                kmer_spec.k
            );
            let motif = decode_ref_kmer(kmer, kmer_spec);
            // N-containing motifs are outside the public A/C/G/T reference k-mer axis
            if motif.contains('N') {
                continue;
            }
            let motif = if canonical {
                make_canonical(motif, true, true)
            } else {
                motif
            };
            // Canonicalization can make multiple encoded keys land in the same output label
            *label_bins[row_idx].entry(motif).or_insert(0.0) += value;
        }
    }

    let motif_order = if all_motifs {
        build_all_ref_kmer_order(kmer_spec, canonical)?
    } else {
        // Observed-only output still uses a sorted axis so repeated runs are stable
        collect_ref_kmer_order(&label_bins)
    };

    // Translate labels to final column indices once, then use integer keys in the writer path
    let motif_columns: FxHashMap<&str, u32> = motif_order
        .iter()
        .enumerate()
        .map(|(column_idx, motif)| {
            Ok((
                motif.as_str(),
                u32::try_from(column_idx)
                    .context("reference k-mer output column index does not fit in u32")?,
            ))
        })
        .collect::<Result<_>>()?;

    // Replace string-keyed row maps with column-index-keyed maps expected by dense and sparse writers
    let count_bins = label_bins
        .into_iter()
        .map(|bin| {
            let mut indexed_bin = FxHashMap::default();
            for (motif, count) in bin {
                let column_idx = motif_columns
                    .get(motif.as_str())
                    .copied()
                    .with_context(|| {
                        format!("missing output column for reference k-mer '{motif}'")
                    })?;
                indexed_bin.insert(column_idx, count);
            }
            Ok(indexed_bin)
        })
        .collect::<Result<Vec<_>>>()?;

    // Normalize at the final column-indexed stage so both storage modes share the same semantics
    let frequency_bins = normalize_count_bins_to_frequencies(count_bins)?;

    Ok((frequency_bins, motif_order))
}

/// Convert count bins to frequencies while retaining one row scaling factor per row.
///
/// The public matrix stores frequencies, not raw counts. For each row, this function sums the
/// retained counts, divides each retained count by that row total, and stores the total separately
/// as `row_scaling_factor`. Downstream code can reconstruct counts with
/// `frequency * row_scaling_factor[row]`.
///
/// Empty rows are valid. They get an empty sparse row and a scaling factor of `0.0`, which avoids a
/// fake denominator while keeping row count and metadata arrays aligned.
///
/// Parameters
/// ----------
/// - `count_bins`:
///   Sparse per-row counts keyed by final output column index
///
/// Returns
/// -------
/// - `Result<RefKmerFrequencyBins>`:
///   Sparse per-row frequencies plus one scaling factor per row
pub(crate) fn normalize_count_bins_to_frequencies(
    count_bins: Vec<RefKmerCountBin>,
) -> Result<RefKmerFrequencyBins> {
    let mut frequency_bins = Vec::with_capacity(count_bins.len());
    let mut row_scaling_factors = Vec::with_capacity(count_bins.len());

    for count_bin in count_bins {
        // Validate and keep counts before computing the denominator used for this row
        let mut kept_counts = Vec::with_capacity(count_bin.len());
        let mut row_total = 0.0;
        for (target_idx, count) in count_bin {
            if !KmerCounts::should_store_weight(count)? {
                continue;
            }
            row_total += count;
            kept_counts.push((target_idx, count));
        }

        if kept_counts.is_empty() {
            // A row with no retained A/C/G/T motifs has no denominator to preserve
            row_scaling_factors.push(0.0);
            frequency_bins.push(FxHashMap::default());
            continue;
        }
        ensure!(
            row_total.is_finite() && row_total > 0.0,
            "reference k-mer row scaling factor {row_total} is not a positive finite value"
        );

        let mut frequency_bin = FxHashMap::default();
        for (target_idx, count) in kept_counts {
            let frequency = count / row_total;
            ensure!(
                frequency.is_finite() && frequency >= 0.0,
                "reference k-mer frequency {frequency} is not a non-negative finite value"
            );
            if frequency > 0.0 {
                // Sparse storage omits exact zeroes but preserves every positive contribution
                frequency_bin.insert(target_idx, frequency);
            }
        }
        row_scaling_factors.push(row_total);
        frequency_bins.push(frequency_bin);
    }

    Ok(RefKmerFrequencyBins {
        frequency_bins,
        row_scaling_factors,
    })
}

/// Return whether the final motif axis already has every expected k-mer label.
///
/// Concrete motif labels are validated by package validation before writing. This helper only
/// answers whether a valid concrete motif axis has the mathematical size of the non-canonical or
/// canonical k-mer set, without generating every motif just to decide the storage mode.
///
/// Grouped motifs from a motifs file can never be treated as a complete k-mer axis because a
/// column may represent more than one motif.
pub(crate) fn ref_kmer_axis_is_complete(
    kmer_size: u8,
    canonical: bool,
    column_kind: SelectedMotifColumnKind,
    motif_labels: &[String],
) -> bool {
    if column_kind != SelectedMotifColumnKind::Motif {
        // Motif-group columns do not have a one-label-per-k-mer interpretation
        return false;
    }
    let Some(expected_num_labels) = complete_ref_kmer_axis_len(usize::from(kmer_size), canonical)
    else {
        // Invalid or oversized all-motif axes should stay on the sparse path
        return false;
    };
    if motif_labels.len() != expected_num_labels {
        // The length check avoids scanning labels when the axis cannot be complete
        return false;
    }

    true
}

/// Return the number of labels in the complete reference k-mer axis.
///
/// Non-canonical output has `4^k` labels. Canonical output groups each k-mer with its reverse
/// complement. Most groups contain two k-mers, but self reverse-complements such as `AT` and `CG`
/// are fixed points and already form a complete group by themselves.
fn complete_ref_kmer_axis_len(kmer_size: usize, canonical: bool) -> Option<usize> {
    if kmer_size == 0 {
        return None;
    }

    let total_kmers = count_possible_acgt_kmers(kmer_size)?;
    if !canonical {
        return total_kmers.try_into().ok();
    }

    let self_reverse_complement_kmers = count_self_reverse_complement_kmers(kmer_size)?;
    let paired_kmers = total_kmers.checked_sub(self_reverse_complement_kmers)?;
    let reverse_complement_pairs = paired_kmers.checked_div(2)?;
    let canonical_kmers = reverse_complement_pairs.checked_add(self_reverse_complement_kmers)?;

    canonical_kmers.try_into().ok()
}

/// Count all possible concrete A/C/G/T k-mers of one length.
fn count_possible_acgt_kmers(kmer_size: usize) -> Option<u64> {
    let exponent: u32 = kmer_size.try_into().ok()?;
    4_u64.checked_pow(exponent)
}

/// Count k-mers that are equal to their own reverse complement.
///
/// Odd k has no fixed points because the middle base would need to complement itself. For even k,
/// choosing the left half determines the right half. For example, choosing `AC` as the left half
/// of a 4-mer forces the right half to be `GT`, so `ACGT` is a self reverse-complement. That leaves
/// `k / 2` free positions, with four base choices at each position, so there are `4^(k / 2)` fixed
/// points.
fn count_self_reverse_complement_kmers(kmer_size: usize) -> Option<u64> {
    if kmer_size % 2 != 0 {
        return Some(0);
    }

    count_possible_acgt_kmers(kmer_size / 2)
}

/// Collect observed motif labels in deterministic sorted order.
///
/// This is the observed-only axis path. The input row maps use hash maps for fast aggregation, so
/// a `BTreeSet` is used here to make the public motif order independent of hash iteration order and
/// tile merge order.
pub(crate) fn collect_ref_kmer_order(bins: &[FxHashMap<String, f64>]) -> Vec<String> {
    let mut motifs = BTreeSet::new();
    for bin in bins {
        motifs.extend(bin.keys().cloned());
    }
    motifs.into_iter().collect()
}

/// Build the full dense reference k-mer axis for `--all-motifs`.
///
/// The full motif set is generated from the same `KmerSpec` used to decode observed counts so label
/// spelling stays consistent across observed-only and all-motif output. Canonical output keeps only
/// unique reverse-complement representatives and does not allocate the pre-collapse motif vector.
pub(crate) fn build_all_ref_kmer_order(
    kmer_spec: &KmerSpec,
    canonical: bool,
) -> Result<Vec<String>> {
    let mut specs = FxHashMap::default();
    let kmer_size: u8 = kmer_spec
        .k
        .try_into()
        .context("k-mer size does not fit in u8 for motif enumeration")?;
    specs.insert(kmer_size, kmer_spec.clone());

    if !canonical {
        return Ok(all_ref_kmer_motifs(kmer_spec.k, &specs));
    }

    let motif_count = 4_u64
        .checked_pow(kmer_spec.k as u32)
        .context("all reference k-mer set overflows u64")?;
    let mut canonical_motifs = BTreeSet::new();
    for radix4_code in 0..motif_count {
        let motif = kmer_spec.decode_kmer(acgt_radix5_code_from_radix4(radix4_code, kmer_spec.k));
        canonical_motifs.insert(make_canonical(motif, true, true));
    }
    Ok(canonical_motifs.into_iter().collect())
}

/// Guard all-motif dense output before counting starts.
///
/// This uses the final public motif axis size. Canonical output therefore checks the number of
/// reverse-complement classes, not the larger pre-canonical `4^k` motif count.
pub(crate) fn ensure_all_ref_kmers_output_size(
    kmer_size: usize,
    canonical: bool,
    n_windows: usize,
) -> Result<()> {
    let n_motifs = complete_ref_kmer_axis_len(kmer_size, canonical)
        .context("all reference k-mer set does not fit in usize")?;

    ensure_dense_ref_kmer_output_size(n_windows, n_motifs).with_context(|| {
        format!("refusing to output all reference k-mers for kmer_size={kmer_size}")
    })
}

/// Guard dense output by the actual matrix size in bytes.
///
/// Dense reference k-mer output stores one `f64` per row and motif. Sparse output is allowed to be
/// much larger on disk if the user really observed many values, but dense output is intentionally
/// capped because it allocates and writes the entire rectangular matrix, including zeroes.
pub(crate) fn ensure_dense_ref_kmer_output_size(n_windows: usize, n_motifs: usize) -> Result<()> {
    let max_dense_output_bytes = max_dense_ref_kmer_output_bytes()?;
    let n_values = (n_windows as u64)
        .checked_mul(n_motifs as u64)
        .context("dense reference k-mer output shape overflows u64")?;
    let bytes = n_values
        .checked_mul(std::mem::size_of::<f64>() as u64)
        .context("dense reference k-mer output byte size overflows u64")?;

    ensure!(
        bytes <= max_dense_output_bytes,
        "Dense reference k-mer output would require {:.2} GiB for {} windows x {} motifs. \
         Reduce the motif space or window count, or set {} to a larger byte limit.",
        bytes as f64 / (1024.0 * 1024.0 * 1024.0),
        n_windows,
        n_motifs,
        MAX_DENSE_REF_KMER_OUTPUT_BYTES_ENV
    );

    Ok(())
}

/// Build grouped row metadata in count-row order.
pub(crate) fn grouped_ref_kmer_row_metadata<'a>(
    group_idx_to_name: &'a FxHashMap<u64, String>,
    chromosomes: &[String],
    grouped_windows_map: &FxHashMap<String, GroupedWindows>,
    blacklist_map: &FxHashMap<String, Vec<Interval<u64>>>,
) -> Result<Vec<RefKmerGroupSummary<'a>>> {
    let mut eligible_windows_by_group: FxHashMap<u64, usize> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0usize))
        .collect();
    let mut total_bp_by_group: FxHashMap<u64, u64> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0u64))
        .collect();
    let mut blacklisted_bp_by_group: FxHashMap<u64, f64> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0.0f64))
        .collect();

    for chromosome in chromosomes {
        let windows = grouped_windows_map
            .get(chromosome)
            .map(|windows| windows.windows_as_slice())
            .unwrap_or(&[]);
        let blacklist_intervals = blacklist_map
            .get(chromosome)
            .map(|intervals| intervals.as_slice())
            .unwrap_or(&[]);
        let mut blacklist_ptr = 0usize;
        for window in windows {
            let (start, end, group_idx) = window.as_tuple();
            ensure!(
                group_idx_to_name.contains_key(&group_idx),
                "grouped reference k-mer window references group_idx {group_idx}, but no group name exists for it"
            );
            let window_bp = end
                .checked_sub(start)
                .context("grouped reference k-mer window end must be >= start")?;
            let blacklist_fraction = compute_blacklist_overlap(
                blacklist_intervals,
                Interval::new(start, end)?,
                0,
                &mut blacklist_ptr,
            );
            *eligible_windows_by_group
                .get_mut(&group_idx)
                .expect("group_idx was validated above") += 1;
            *total_bp_by_group
                .get_mut(&group_idx)
                .expect("group_idx was validated above") += window_bp;
            *blacklisted_bp_by_group
                .get_mut(&group_idx)
                .expect("group_idx was validated above") += blacklist_fraction * window_bp as f64;
        }
    }

    let mut summaries = Vec::with_capacity(group_idx_to_name.len());
    for row_index in 0..group_idx_to_name.len() {
        let group_idx = row_index as u64;
        let Some(group_name) = group_idx_to_name.get(&group_idx) else {
            let mut observed_indices: Vec<u64> = group_idx_to_name.keys().copied().collect();
            observed_indices.sort_unstable();
            bail!(
                "reference k-mer group indices must match count rows 0..{} but observed {:?}",
                group_idx_to_name.len().saturating_sub(1),
                observed_indices
            );
        };
        validate_zarr_label(group_name, "group_name")?;
        let total_bp = total_bp_by_group.get(&group_idx).copied().unwrap_or(0);
        let blacklisted_fraction = if total_bp == 0 {
            0.0
        } else {
            blacklisted_bp_by_group
                .get(&group_idx)
                .copied()
                .unwrap_or(0.0)
                / total_bp as f64
        };
        summaries.push(RefKmerGroupSummary {
            group_idx,
            group_name,
            eligible_windows: eligible_windows_by_group
                .get(&group_idx)
                .copied()
                .unwrap_or(0),
            blacklisted_fraction,
        });
    }

    Ok(summaries)
}

/// Write a complete reference k-mer frequency package as a Zarr V3 store.
///
/// The sparse frequency rows are expected to be the normalized output of
/// `normalize_count_bins_to_frequencies` or an equivalent upstream step.
pub(crate) fn write_ref_kmer_zarr(
    store_path: &Path,
    package: RefKmerZarrPackage<'_>,
) -> Result<()> {
    validate_ref_kmer_package(&package)?;
    let store = create_zarr_store(store_path, "reference k-mer")?;
    let storage_mode = if package.write_dense_output {
        "dense"
    } else {
        "sparse_coo"
    };
    let row_mode = row_mode_name(&package.row_metadata);
    write_root_metadata(store.clone(), storage_mode, row_mode, &package)?;
    write_motif_metadata(store.clone(), package.motif_labels, package.column_kind)?;
    write_row_metadata(
        store.clone(),
        package.row_metadata,
        package.frequency_bins.len(),
    )?;
    write_row_scaling_factors(store.clone(), package.row_scaling_factors)?;

    if package.write_dense_output {
        let frequencies =
            stack_ref_kmer_frequencies(package.frequency_bins, package.motif_labels.len())?;
        write_dense_frequencies(store.clone(), frequencies.view())?;
    } else {
        write_sparse_frequencies(
            store.clone(),
            package.frequency_bins,
            package.motif_labels.len(),
        )?;
    }

    let reference_contig_footprint_json = serde_json::to_vec(package.reference_contig_footprint)?;
    write_single_chunk_zarr_array(
        store,
        "reference_contig_footprint_json",
        &[reference_contig_footprint_json.len()],
        &["json_byte"],
        &reference_contig_footprint_json,
        data_type::uint8(),
        0u8,
        json!({"long_name": "JSON-encoded reference contig footprint"}),
    )?;

    Ok(())
}

/// Convenience path builder for callers that write into a temporary output directory.
pub(crate) fn ref_kmer_zarr_path(output_dir: &Path, prefix: &str) -> PathBuf {
    output_dir.join(dot_join(&[prefix, "ref_kmer_counts.zarr"]))
}

fn decode_ref_kmer(kmer: Kmer, kmer_spec: &KmerSpec) -> String {
    let motif = kmer_spec.decode_kmer(kmer.code);
    match kmer.orientation {
        KmerOrientation::Forward => motif,
        KmerOrientation::Reverse => crate::shared::base::rev_complement(&motif),
    }
}

fn max_dense_ref_kmer_output_bytes() -> Result<u64> {
    match env::var(MAX_DENSE_REF_KMER_OUTPUT_BYTES_ENV) {
        Ok(raw_value) => {
            let parsed = raw_value.parse::<u64>().with_context(|| {
                format!(
                    "{} must be a positive integer byte count, got '{}'",
                    MAX_DENSE_REF_KMER_OUTPUT_BYTES_ENV, raw_value
                )
            })?;
            ensure!(
                parsed > 0,
                "{} must be > 0 bytes",
                MAX_DENSE_REF_KMER_OUTPUT_BYTES_ENV
            );
            Ok(parsed)
        }
        Err(env::VarError::NotPresent) => Ok(DEFAULT_MAX_DENSE_REF_KMER_OUTPUT_BYTES),
        Err(env::VarError::NotUnicode(_)) => bail!(
            "{} must contain valid UTF-8 text",
            MAX_DENSE_REF_KMER_OUTPUT_BYTES_ENV
        ),
    }
}

fn validate_ref_kmer_package(package: &RefKmerZarrPackage<'_>) -> Result<()> {
    ensure!(
        !package.frequency_bins.is_empty(),
        "reference k-mer Zarr output requires at least one row"
    );
    validate_ref_kmer_motif_axis(
        package.kmer_size,
        package.canonical,
        package.column_kind,
        package.motif_labels,
    )?;
    ensure!(
        package.row_scaling_factors.len() == package.frequency_bins.len(),
        "reference k-mer row scaling factors ({}) did not match frequency rows ({})",
        package.row_scaling_factors.len(),
        package.frequency_bins.len()
    );
    for (row_idx, scaling_factor) in package.row_scaling_factors.iter().enumerate() {
        ensure!(
            scaling_factor.is_finite() && *scaling_factor >= 0.0,
            "reference k-mer row scaling factor at row {} is not a non-negative finite value: {}",
            row_idx,
            scaling_factor
        );
    }
    if package.write_dense_output {
        ensure!(
            !package.motif_labels.is_empty(),
            "dense reference k-mer output requires at least one motif column"
        );
        ensure_dense_ref_kmer_output_size(
            package.frequency_bins.len(),
            package.motif_labels.len(),
        )?;
    }
    Ok(())
}

fn validate_ref_kmer_motif_axis(
    kmer_size: u8,
    canonical: bool,
    column_kind: SelectedMotifColumnKind,
    motif_labels: &[String],
) -> Result<()> {
    ensure!(
        kmer_size > 0,
        "reference k-mer Zarr output requires positive k-mer size"
    );
    match column_kind {
        SelectedMotifColumnKind::Motif => {
            validate_concrete_ref_kmer_motif_axis(kmer_size, canonical, motif_labels)
        }
        SelectedMotifColumnKind::MotifGroup => validate_ref_kmer_motif_group_axis(motif_labels),
    }
}

fn validate_concrete_ref_kmer_motif_axis(
    kmer_size: u8,
    canonical: bool,
    motif_labels: &[String],
) -> Result<()> {
    let expected_len = usize::from(kmer_size);
    let mut seen = BTreeSet::new();
    for (motif_index, motif) in motif_labels.iter().enumerate() {
        validate_zarr_label(motif, "motif")?;
        ensure!(
            motif.len() == expected_len,
            "reference k-mer motif label `{}` at motif index {} has length {}, expected {}",
            motif,
            motif_index,
            motif.len(),
            expected_len
        );
        for base in motif.bytes() {
            ensure!(
                matches!(base, b'A' | b'C' | b'G' | b'T'),
                "reference k-mer motif label `{}` at motif index {} contains invalid base `{}`",
                motif,
                motif_index,
                base as char
            );
        }
        if canonical {
            let canonical_motif = make_canonical(motif.clone(), true, true);
            ensure!(
                canonical_motif == *motif,
                "canonical reference k-mer motif label `{}` at motif index {} should be represented as `{}`",
                motif,
                motif_index,
                canonical_motif
            );
        }
        ensure!(
            seen.insert(motif.as_str()),
            "duplicate reference k-mer motif label `{}` at motif index {}",
            motif,
            motif_index
        );
    }
    Ok(())
}

fn validate_ref_kmer_motif_group_axis(motif_groups: &[String]) -> Result<()> {
    let mut seen = BTreeSet::new();
    for (motif_index, motif_group) in motif_groups.iter().enumerate() {
        validate_zarr_label(motif_group, "motif_group")?;
        ensure!(
            seen.insert(motif_group.as_str()),
            "duplicate reference k-mer motif-group label `{}` at motif index {}",
            motif_group,
            motif_index
        );
    }
    Ok(())
}

fn row_mode_name(row_metadata: &RefKmerRowMetadata<'_>) -> &'static str {
    match row_metadata {
        RefKmerRowMetadata::Global => "global",
        RefKmerRowMetadata::Windows { row_mode, .. } => match row_mode {
            RefKmerWindowRowMode::Size => "size",
            RefKmerWindowRowMode::Bed => "bed",
        },
        RefKmerRowMetadata::Groups(_) => "grouped_bed",
    }
}

fn write_root_metadata(
    store: Arc<FilesystemStore>,
    storage_mode: &str,
    row_mode: &str,
    package: &RefKmerZarrPackage<'_>,
) -> Result<()> {
    let mut attributes = json!({
        "cfdnalab_schema": "ref_kmer_frequencies",
        "cfdnalab_schema_version": CFDNALAB_REF_KMER_SCHEMA_VERSION,
        "storage_mode": storage_mode,
        "row_mode": row_mode,
        "motif_axis_kind": motif_column_kind_name(package.column_kind),
        "value_units": "reference_kmer_frequency",
        "count_units": "reference_kmer_count",
        "row_scaling_factor_array": "row_scaling_factor",
        "count_reconstruction": "reference_kmer_count = frequency * row_scaling_factor[row]",
        "kmer_size": package.kmer_size,
        "canonical": package.canonical,
        "all_motifs": package.all_motifs,
        "assign_by": assign_by_name(package.assign_by),
        "primary_array": null,
        "primary_group": null,
    });
    if storage_mode == "dense" {
        attributes["primary_array"] = json!("frequencies");
    } else if storage_mode == "sparse_coo" {
        attributes["primary_group"] = json!("sparse");
        attributes["sparse_format"] = json!("coo");
        attributes["sparse_indices_base"] = json!(0);
    }
    write_zarr_root_metadata(store, "reference k-mer", attributes)
}

fn motif_column_kind_name(column_kind: SelectedMotifColumnKind) -> &'static str {
    match column_kind {
        SelectedMotifColumnKind::Motif => "motif",
        SelectedMotifColumnKind::MotifGroup => "motif_group",
    }
}

fn assign_by_name(assign_by: WindowAssigner) -> String {
    match assign_by {
        WindowAssigner::CountOverlap => "count-overlap".to_string(),
        WindowAssigner::Any => "any".to_string(),
        WindowAssigner::All => "all".to_string(),
        WindowAssigner::Midpoint => "midpoint".to_string(),
        WindowAssigner::Proportion(threshold) => format!("proportion={threshold}"),
    }
}

fn write_motif_metadata(
    store: Arc<FilesystemStore>,
    labels: &[String],
    column_kind: SelectedMotifColumnKind,
) -> Result<()> {
    match column_kind {
        SelectedMotifColumnKind::Motif => write_motif_label_metadata(store, labels),
        SelectedMotifColumnKind::MotifGroup => write_motif_group_metadata(store, labels),
    }
}

fn write_motif_label_metadata(store: Arc<FilesystemStore>, motifs: &[String]) -> Result<()> {
    let motif_width = validated_motif_width(motifs)?;
    let motif_axis = checked_index_axis(motifs.len(), "motif")?;
    let motif_byte_axis = checked_index_axis(motif_width, "motif_byte")?;
    let motif_ascii = encode_motif_ascii(motifs, motif_width);

    write_single_chunk_zarr_array(
        store.clone(),
        "motif_index",
        &[motifs.len()],
        &["motif"],
        &motif_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based motif column index",
            "label_array": "motif_ascii",
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "motif_byte",
        &[motif_width],
        &["motif_byte"],
        &motif_byte_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "zero-based byte offset within reference k-mer labels" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "motif_ascii",
        &[motifs.len(), motif_width],
        &["motif", "motif_byte"],
        &motif_ascii,
        data_type::uint8(),
        ZARR_ASCII_FILL_VALUE,
        json!({
            "long_name": "fixed-width ASCII reference k-mer labels",
            "description": "Decode each [motif, motif_byte] row as ASCII to recover the motif label.",
        }),
    )?;
    Ok(())
}

fn write_motif_group_metadata(store: Arc<FilesystemStore>, motif_groups: &[String]) -> Result<()> {
    for motif_group in motif_groups {
        validate_zarr_label(motif_group, "motif_group")?;
    }
    let motif_axis = checked_index_axis(motif_groups.len(), "motif")?;
    let motif_group_refs: Vec<&str> = motif_groups.iter().map(String::as_str).collect();

    write_single_chunk_zarr_array(
        store,
        "motif_index",
        &[motif_groups.len()],
        &["motif"],
        &motif_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based motif group column index",
            "label_field": "motif_group",
            "labels": motif_group_refs,
        }),
    )?;
    Ok(())
}

fn validated_motif_width(motifs: &[String]) -> Result<usize> {
    let Some(first_motif) = motifs.first() else {
        return Ok(0);
    };
    validate_motif_ascii(first_motif)?;
    let motif_width = first_motif.len();
    for motif in &motifs[1..] {
        validate_motif_ascii(motif)?;
        ensure!(
            motif.len() == motif_width,
            "reference k-mer labels must have one fixed ASCII width, but '{}' has {} bytes and expected {}",
            motif,
            motif.len(),
            motif_width
        );
    }
    Ok(motif_width)
}

fn validate_motif_ascii(motif: &str) -> Result<()> {
    validate_zarr_label(motif, "motif")?;
    ensure!(
        motif.is_ascii(),
        "reference k-mer labels must be ASCII to be stored in motif_ascii"
    );
    Ok(())
}

fn encode_motif_ascii(motifs: &[String], motif_width: usize) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(motifs.len() * motif_width);
    for motif in motifs {
        bytes.extend_from_slice(motif.as_bytes());
    }
    bytes
}

fn write_row_metadata(
    store: Arc<FilesystemStore>,
    row_metadata: RefKmerRowMetadata<'_>,
    n_rows: usize,
) -> Result<()> {
    ensure!(
        n_rows > 0,
        "reference k-mer Zarr output requires at least one row"
    );
    let row = checked_index_axis(n_rows, "row")?;
    let mut row_attributes = json!({ "long_name": "zero-based frequency row index" });
    if matches!(&row_metadata, RefKmerRowMetadata::Global) {
        row_attributes["label_field"] = json!("row_label");
        row_attributes["labels"] = json!(["global"]);
    }
    write_single_chunk_zarr_array(
        store.clone(),
        "row",
        &[n_rows],
        &["row"],
        &row,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        row_attributes,
    )?;

    match row_metadata {
        RefKmerRowMetadata::Global => {
            ensure!(
                n_rows == 1,
                "global reference k-mer output must have one row"
            );
        }
        RefKmerRowMetadata::Windows { bin_info, .. } => {
            write_window_row_metadata(store, bin_info, n_rows)?;
        }
        RefKmerRowMetadata::Groups(groups) => {
            write_group_row_metadata(store, &groups, n_rows)?;
        }
    }
    Ok(())
}

fn write_window_row_metadata(
    store: Arc<FilesystemStore>,
    bin_info: &[WindowBinInfo],
    n_rows: usize,
) -> Result<()> {
    ensure!(
        bin_info.len() == n_rows,
        "reference k-mer window metadata rows ({}) did not match frequency rows ({})",
        bin_info.len(),
        n_rows
    );
    for pair in bin_info.windows(2) {
        ensure!(
            pair[0].output_index < pair[1].output_index,
            "reference k-mer window metadata must be sorted by increasing output_index"
        );
    }

    let chromosome_names = chromosome_axis_from_windows(bin_info)?;
    let chromosome_name_refs: Vec<&str> = chromosome_names.iter().map(String::as_str).collect();
    let chromosome_axis = checked_index_axis(chromosome_names.len(), "chromosome")?;
    let chromosome_lookup: FxHashMap<&str, usize> = chromosome_names
        .iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect();
    let row_chromosome = bin_info
        .iter()
        .map(|entry| {
            let index = chromosome_lookup
                .get(entry.chromosome.as_str())
                .copied()
                .context("window chromosome missing from chromosome axis")?;
            checked_i32(index, "row_chromosome")
        })
        .collect::<Result<Vec<_>>>()?;
    let row_start_bp: Vec<i64> = bin_info
        .iter()
        .map(|entry| checked_i64(entry.start, "row_start_bp"))
        .collect::<Result<_>>()?;
    let row_end_bp: Vec<i64> = bin_info
        .iter()
        .map(|entry| checked_i64(entry.end, "row_end_bp"))
        .collect::<Result<_>>()?;
    let blacklisted_fraction: Vec<f64> = bin_info
        .iter()
        .map(|entry| entry.blacklisted_fraction)
        .collect();

    write_single_chunk_zarr_array(
        store.clone(),
        "chromosome",
        &[chromosome_names.len()],
        &["chromosome"],
        &chromosome_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based chromosome index",
            "label_field": "chromosome_name",
            "labels": chromosome_name_refs,
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "row_chromosome",
        &[n_rows],
        &["row"],
        &row_chromosome,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "chromosome index for each frequency row" }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "row_start_bp",
        &[n_rows],
        &["row"],
        &row_start_bp,
        data_type::int64(),
        ZARR_INT64_FILL_VALUE,
        json!({ "long_name": "inclusive row start coordinate", "units": "bp" }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "row_end_bp",
        &[n_rows],
        &["row"],
        &row_end_bp,
        data_type::int64(),
        ZARR_INT64_FILL_VALUE,
        json!({ "long_name": "exclusive row end coordinate", "units": "bp" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "blacklisted_fraction",
        &[n_rows],
        &["row"],
        &blacklisted_fraction,
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({ "long_name": "fraction of row coordinates overlapping blacklist intervals" }),
    )?;

    Ok(())
}

fn write_group_row_metadata(
    store: Arc<FilesystemStore>,
    groups: &[RefKmerGroupSummary<'_>],
    n_rows: usize,
) -> Result<()> {
    ensure!(
        groups.len() == n_rows,
        "reference k-mer group metadata rows ({}) did not match frequency rows ({})",
        groups.len(),
        n_rows
    );
    for group in groups {
        validate_zarr_label(group.group_name, "group_name")?;
    }
    let group: Vec<i32> = groups
        .iter()
        .map(|group| checked_i32(group.group_idx, "group_idx"))
        .collect::<Result<Vec<_>>>()?;
    let group_name_refs: Vec<&str> = groups.iter().map(|group| group.group_name).collect();
    let eligible_windows: Vec<i32> = groups
        .iter()
        .map(|group| checked_i32(group.eligible_windows, "eligible_windows"))
        .collect::<Result<Vec<_>>>()?;
    let blacklisted_fraction: Vec<f64> = groups
        .iter()
        .map(|group| group.blacklisted_fraction)
        .collect();

    write_single_chunk_zarr_array(
        store.clone(),
        "group",
        &[n_rows],
        &["row"],
        &group,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based group index",
            "description": "For grouped BED output, this matches the frequency row index.",
            "label_field": "group_name",
            "labels": group_name_refs,
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "eligible_windows",
        &[n_rows],
        &["row"],
        &eligible_windows,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "eligible grouped BED windows retained per group" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "blacklisted_fraction",
        &[n_rows],
        &["row"],
        &blacklisted_fraction,
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({ "long_name": "length-weighted blacklisted fraction across group windows" }),
    )?;

    Ok(())
}

fn write_row_scaling_factors(
    store: Arc<FilesystemStore>,
    row_scaling_factors: &[f64],
) -> Result<()> {
    write_single_chunk_zarr_array(
        store,
        "row_scaling_factor",
        &[row_scaling_factors.len()],
        &["row"],
        row_scaling_factors,
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({
            "long_name": "reference k-mer counts per row",
            "units": "reference_kmer_count",
            "description": "Multiply row frequencies by this value to reconstruct counts.",
        }),
    )
}

fn stack_ref_kmer_frequencies(
    bins: &[RefKmerFrequencyBin],
    n_columns: usize,
) -> Result<Array2<f64>> {
    let mut frequencies = Array2::<f64>::zeros((bins.len(), n_columns));

    for (row, bin) in bins.iter().enumerate() {
        for (&target_idx, &frequency) in bin {
            let column = target_idx as usize;
            ensure!(
                column < n_columns,
                "reference k-mer column index {} is out of bounds for {} columns",
                column,
                n_columns
            );
            frequencies[(row, column)] = frequency;
        }
    }

    Ok(frequencies)
}

fn write_dense_frequencies(
    store: Arc<FilesystemStore>,
    frequencies: ArrayView2<'_, f64>,
) -> Result<()> {
    ensure!(
        frequencies.shape().iter().all(|dimension| *dimension > 0),
        "dense reference k-mer frequencies cannot have an empty axis"
    );
    let shape = [frequencies.shape()[0], frequencies.shape()[1]];
    let chunk_shape = dense_frequency_chunk_shape(shape)?;
    let array = create_zarr_array(
        store,
        "frequencies",
        &shape,
        &chunk_shape,
        &["row", "motif"],
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({
            "long_name": "reference k-mer frequency",
            "units": "reference_kmer_frequency",
        }),
    )?;
    let owned_frequencies = frequencies.as_standard_layout();
    let values = owned_frequencies
        .as_slice()
        .context("standard-layout reference k-mer frequencies were not contiguous")?;
    array
        .store_array_subset(&array.subset_all(), values)
        .context("write dense reference k-mer Zarr frequencies")?;
    Ok(())
}

fn dense_frequency_chunk_shape(shape: [usize; 2]) -> Result<[usize; 2]> {
    let total_cells = shape
        .iter()
        .try_fold(1usize, |cells, dimension| cells.checked_mul(*dimension))
        .context("reference k-mer dense Zarr frequency shape overflow")?;
    if total_cells <= TARGET_DENSE_FREQUENCY_CHUNK_CELLS {
        return Ok(shape);
    }

    let motif_chunk = shape[1].clamp(1, TARGET_DENSE_FREQUENCY_CHUNK_CELLS);
    let row_chunk = shape[0]
        .min((TARGET_DENSE_FREQUENCY_CHUNK_CELLS / motif_chunk).max(1))
        .max(1);
    Ok([row_chunk, motif_chunk])
}

/// Write normalized reference k-mer frequencies in sparse COO form.
///
/// The sparse package stores parallel `row`, `motif`, and `frequency` vectors plus the dense
/// matrix shape they represent. This function uses a count pass followed by a chunked write pass
/// so array metadata has the exact `nnz` length without collecting every COO entry in memory.
/// Frequency values are expected to have been normalized before this writer is called.
fn write_sparse_frequencies(
    store: Arc<FilesystemStore>,
    bins: &[RefKmerFrequencyBin],
    n_columns: usize,
) -> Result<()> {
    // Count first so validation fails before sparse arrays are created
    let nnz = count_sparse_entries(bins, n_columns)?;
    write_zarr_group_metadata(
        store.clone(),
        "/sparse",
        "reference k-mer sparse frequencies",
        json!({
            "long_name": "Sparse COO reference k-mer frequency arrays",
            "sparse_format": "coo",
            "sparse_indices_base": 0,
        }),
    )?;

    let shape = vec![
        checked_i32(bins.len(), "sparse row count")?,
        checked_i32(n_columns, "sparse motif count")?,
    ];
    let sparse_dimension = checked_index_axis(2, "sparse_dimension")?;
    let vector_shape = [nnz];
    let vector_chunk_shape = [sparse_coo_chunk_entries(nnz)];

    let row_array = create_zarr_array(
        store.clone(),
        "sparse/row",
        &vector_shape,
        &vector_chunk_shape,
        &["nnz"],
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "COO row index" }),
    )?;
    let motif_array = create_zarr_array(
        store.clone(),
        "sparse/motif",
        &vector_shape,
        &vector_chunk_shape,
        &["nnz"],
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "COO motif index" }),
    )?;
    let frequency_array = create_zarr_array(
        store.clone(),
        "sparse/frequency",
        &vector_shape,
        &vector_chunk_shape,
        &["nnz"],
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({
            "long_name": "reference k-mer frequency",
            "units": "reference_kmer_frequency",
        }),
    )?;
    write_sparse_entry_chunks(
        &row_array,
        &motif_array,
        &frequency_array,
        bins,
        n_columns,
        nnz,
    )?;

    write_single_chunk_zarr_array(
        store.clone(),
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &shape,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "dense shape represented by sparse COO arrays" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "sparse/sparse_dimension",
        &[2],
        &["sparse_dimension"],
        &sparse_dimension,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based sparse shape dimension index",
            "label_field": "sparse_dimension_name",
            "labels": ["row", "motif"],
        }),
    )?;
    Ok(())
}

/// Count non-zero entries for the sparse COO arrays.
///
/// Sparse Zarr arrays need their final one-dimensional shape before they are created. This pass
/// counts exact non-zero frequency entries and validates sparse coordinates while leaving
/// frequency normalization assumptions to the upstream normalization step.
fn count_sparse_entries(bins: &[RefKmerFrequencyBin], n_columns: usize) -> Result<usize> {
    let mut nnz = 0usize;
    for (row_index, bin) in bins.iter().enumerate() {
        checked_i32(row_index, "row index")?;
        for (&target_idx, &frequency) in bin {
            let column = target_idx as usize;
            ensure!(
                column < n_columns,
                "sparse reference k-mer column index {} is out of bounds for {} columns",
                column,
                n_columns
            );
            if frequency != 0.0 {
                checked_i32(target_idx, "motif index")?;
                nnz = nnz
                    .checked_add(1)
                    .context("sparse reference k-mer entry count overflow")?;
            }
        }
    }
    Ok(nnz)
}

/// Choose the one-dimensional chunk length for sparse COO vectors.
///
/// Empty sparse matrices have shape `[0]`, but Zarr chunk shapes still need a positive length.
/// Non-empty matrices use at most `TARGET_SPARSE_COO_CHUNK_ENTRIES` entries per chunk.
fn sparse_coo_chunk_entries(nnz: usize) -> usize {
    nnz.max(1).min(TARGET_SPARSE_COO_CHUNK_ENTRIES)
}

/// Stream sparse rows into aligned COO chunks.
///
/// Output is deterministic row-major COO order: rows follow the input row order, and motif indices
/// are sorted within each row. The writer keeps memory bounded by one output chunk plus the keys
/// for the current row, instead of staging the full sparse matrix.
fn write_sparse_entry_chunks(
    row_array: &Array<FilesystemStore>,
    motif_array: &Array<FilesystemStore>,
    frequency_array: &Array<FilesystemStore>,
    bins: &[RefKmerFrequencyBin],
    n_columns: usize,
    nnz: usize,
) -> Result<()> {
    if nnz == 0 {
        return Ok(());
    }

    let chunk_capacity = sparse_coo_chunk_entries(nnz);
    let mut row_chunk = Vec::with_capacity(chunk_capacity);
    let mut motif_chunk = Vec::with_capacity(chunk_capacity);
    let mut frequency_chunk = Vec::with_capacity(chunk_capacity);
    let mut write_offset = 0usize;

    for (row_index, bin) in bins.iter().enumerate() {
        let row_index = checked_i32(row_index, "row index")?;
        // Sort each row locally to keep COO output deterministic without staging all entries
        let mut motif_indices: Vec<u32> = bin.keys().copied().collect();
        motif_indices.sort_unstable();

        for motif_index in motif_indices {
            let column = motif_index as usize;
            ensure!(
                column < n_columns,
                "sparse reference k-mer column index {} is out of bounds for {} columns",
                column,
                n_columns
            );
            let frequency = bin[&motif_index];
            if frequency == 0.0 {
                continue;
            }

            row_chunk.push(row_index);
            motif_chunk.push(checked_i32(motif_index, "motif index")?);
            frequency_chunk.push(frequency);

            if row_chunk.len() == chunk_capacity {
                write_sparse_entry_chunk(
                    row_array,
                    motif_array,
                    frequency_array,
                    write_offset,
                    &row_chunk,
                    &motif_chunk,
                    &frequency_chunk,
                )?;
                write_offset += row_chunk.len();
                row_chunk.clear();
                motif_chunk.clear();
                frequency_chunk.clear();
            }
        }
    }

    if !row_chunk.is_empty() {
        write_sparse_entry_chunk(
            row_array,
            motif_array,
            frequency_array,
            write_offset,
            &row_chunk,
            &motif_chunk,
            &frequency_chunk,
        )?;
        write_offset += row_chunk.len();
    }
    ensure!(
        write_offset == nnz,
        "wrote {} sparse reference k-mer entries but expected {}",
        write_offset,
        nnz
    );

    Ok(())
}

/// Write one aligned sparse COO slice to the Zarr arrays.
///
/// The three input slices are parallel columns for the same `nnz` range. Each call writes the same
/// subset into `sparse/row`, `sparse/motif`, and `sparse/frequency`, so readers can reconstruct
/// each COO entry by taking the same index from all three arrays.
fn write_sparse_entry_chunk(
    row_array: &Array<FilesystemStore>,
    motif_array: &Array<FilesystemStore>,
    frequency_array: &Array<FilesystemStore>,
    write_offset: usize,
    row_chunk: &[i32],
    motif_chunk: &[i32],
    frequency_chunk: &[f64],
) -> Result<()> {
    ensure!(
        row_chunk.len() == motif_chunk.len() && row_chunk.len() == frequency_chunk.len(),
        "sparse reference k-mer chunk columns have inconsistent lengths"
    );
    let start = u64::try_from(write_offset).context("sparse write offset exceeds u64")?;
    let end_offset = write_offset
        .checked_add(row_chunk.len())
        .context("sparse write offset overflow")?;
    let end = u64::try_from(end_offset).context("sparse write end offset exceeds u64")?;
    let subset = ArraySubset::new_with_ranges(&[start..end]);

    row_array
        .store_array_subset(&subset, row_chunk)
        .context("write sparse reference k-mer Zarr row chunk")?;
    motif_array
        .store_array_subset(&subset, motif_chunk)
        .context("write sparse reference k-mer Zarr motif chunk")?;
    frequency_array
        .store_array_subset(&subset, frequency_chunk)
        .context("write sparse reference k-mer Zarr frequency chunk")?;

    Ok(())
}

fn chromosome_axis_from_windows(bin_info: &[WindowBinInfo]) -> Result<Vec<String>> {
    let mut chromosome_names = Vec::new();
    let mut seen = FxHashMap::default();
    for entry in bin_info {
        validate_zarr_label(&entry.chromosome, "chromosome_name")?;
        if !seen.contains_key(entry.chromosome.as_str()) {
            seen.insert(entry.chromosome.as_str(), chromosome_names.len());
            chromosome_names.push(entry.chromosome.clone());
        }
    }
    ensure!(
        !chromosome_names.is_empty(),
        "windowed reference k-mer output requires at least one chromosome"
    );
    Ok(chromosome_names)
}

#[cfg(test)]
mod tests {
    include!("zarr_tests.rs");
}
