#[cfg(not(feature = "cmd_fragment_kmers"))]
use crate::shared::kmers::kmer_codec::KmerSpec;
#[cfg(feature = "cmd_fragment_kmers")]
use crate::shared::{
    base::make_canonical,
    kmers::kmer_codec::{Kmer, KmerSpec},
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
#[cfg(feature = "cmd_fragment_kmers")]
use fxhash::FxHashSet;

/// Per-k map of “reference” counts
#[derive(Debug, Clone, PartialEq)]
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) struct DecodedCounts {
    pub(crate) counts: FxHashMap<u8, FxHashMap<String, f64>>, // k -> motif -> count
}

/// Prepare decoded counts for all kmer sizes in all windows.
///
/// Extracts motifs per kmer spec to allow future padding.
/// For kmers of size 1..6, this includes all possible motifs.
/// For larger kmer sizes, only the seen motifs is included as the number otherwise explodes.
///
/// * `windows`        – slice of per-window raw counts
/// * `canonical`      – canonical reverse complements when true
/// * `kmer_specs`     – validated specs for every k we want to keep
///
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn prepare_decoded_counts(
    windows: &[DecodedCounts],
    canonical: bool,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
) -> (Vec<DecodedCounts>, FxHashMap<u8, Vec<String>>) {
    let n_windows = windows.len();

    // Initialise one empty DecodedCounts per window
    let mut out = vec![
        DecodedCounts {
            counts: FxHashMap::default()
        };
        n_windows
    ];

    let mut motifs_by_k: FxHashMap<u8, Vec<String>> = FxHashMap::default();

    // Loop over every k we validated
    for &k in kmer_specs.keys() {
        // Reference (match) bins for this k
        let (count_bins, motifs) =
            prepare_kmer_category(windows, kmer_specs, k as usize, canonical, k <= 6);

        // Insert into the corresponding window
        for i in 0..n_windows {
            out[i].counts.insert(k, count_bins[i].clone());
        }
        motifs_by_k.insert(k, motifs);
    }

    (out, motifs_by_k)
}

#[cfg(feature = "cmd_fragment_kmers")]
fn prepare_kmer_category(
    windows: &[DecodedCounts],
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    k: usize,
    canonical: bool,
    ensure_all: bool,
) -> (Vec<FxHashMap<String, f64>>, Vec<String>) {
    // Extract the raw maps
    let raw_bins = extract_bins(windows, k, canonical);

    // Build the (canonical) motif list once, if requested.
    let base_motifs: Vec<String> = if ensure_all {
        all_motifs(k, kmer_specs)
    } else {
        Vec::new()
    };

    // Build the (canonical) motif list *once* so we know what to pad with
    let mut motifs = collect_motifs(&raw_bins, base_motifs, canonical, ensure_all);
    motifs.sort_unstable();

    (raw_bins, motifs)
}

/// Collect per-window bins for the requested motif type and (optionally)
/// canonical them into strand-agnostic form.
///
/// * `windows` – slice of `DecodedCounts` (“one window” each).
/// * `k` – kmer-size to pull out of every `DecodedCounts`.
/// * `canonical` – if `true`, run the appropriate collapse_*_map helper.
///
/// Returns a fresh `Vec<FxHashMap<String, f64>>` – one map per window.
#[cfg(feature = "cmd_fragment_kmers")]
fn extract_bins(
    windows: &[DecodedCounts],
    k: usize, // pattern only; field values are ignored
    canonical: bool,
) -> Vec<FxHashMap<String, f64>> {
    windows
        .iter()
        .map(|dc| {
            // Pick the raw map for this window
            let raw: FxHashMap<String, f64> =
                dc.counts.get(&(k as u8)).cloned().unwrap_or_default();

            // Collapse if requested, otherwise return the raw map
            if canonical { collapse_map(raw) } else { raw }
        })
        .collect()
}

/// Collect motifs for a category, optionally ensuring the complete set and filtering 'N'
#[cfg(feature = "cmd_fragment_kmers")]
fn collect_motifs(
    windows: &[FxHashMap<String, f64>],
    base_motifs: Vec<String>,
    canonical: bool,
    ensure_all: bool,
) -> Vec<String> {
    // Set of motifs to keep
    let set: FxHashSet<String> = if ensure_all {
        base_motifs.into_iter().collect()
    } else {
        windows.iter().flat_map(|m| m.keys().cloned()).collect()
    };

    // Strand-collapse if requested
    let collapsed_set = if canonical { collapse_set(&set) } else { set };

    // Convert to sorted Vec
    let mut v: Vec<String> = collapsed_set.into_iter().collect();
    v.sort_unstable();
    v
}

/// Collapse a map of kmer counts into canonical keys, summing counts
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn collapse_map(map: FxHashMap<String, f64>) -> FxHashMap<String, f64> {
    let mut out: FxHashMap<String, f64> = FxHashMap::default();
    out.reserve(map.len());

    for (kmer, count) in map {
        let canon = make_canonical(kmer, true, false);
        *out.entry(canon).or_insert(0.0) += count;
    }

    out
}

/// Collapse a set of reference kmers into a set of canonical keys
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn collapse_set(set: &FxHashSet<String>) -> FxHashSet<String> {
    set.iter()
        .map(|m| make_canonical(m.to_string(), true, false))
        .collect()
}

/// Return all possible reference motifs (4ᵏ) for a given k.
///
/// No motifs with 'N' are returned.
pub(crate) fn all_motifs(k: usize, specs: &FxHashMap<u8, KmerSpec>) -> Vec<String> {
    let spec = &specs[&(k as u8)];
    let motif_count = 4u64.pow(k as u32);
    (0..motif_count)
        .map(|code| spec.decode_kmer(acgt_radix5_code_from_radix4(code, k)))
        .collect()
}

/// Convert an A/C/G/T-only radix-4 motif index into the shared radix-5 k-mer code.
///
/// The shared [`KmerSpec`] decoder expects codes in the package's ordinary radix-5 encoding:
///
/// - `A = 0`
/// - `C = 1`
/// - `G = 2`
/// - `T = 3`
/// - `N = 4`
///
/// `all_motifs()` only needs the complete A/C/G/T set. A radix-4 counter represents exactly that
/// no-N set, because each digit has only the four allowed A/C/G/T states.
///
/// The function has to keep integer values and positional representations separate. A decimal value
/// is the ordinary base-10 integer value, like the value stored in a `u64`. A digit is one
/// coefficient in a positional representation of that value in a chosen base. The same decimal value
/// can therefore have different digits in different bases. For example, decimal `6` is `12` in
/// radix-4, meaning digits `[1, 2]`, because `1 * 4 + 2 = 6`. The same decimal `6` is `11` in
/// radix-5, meaning digits `[1, 1]`, because `1 * 5 + 1 = 6`.
///
/// The numeric radix-4 value cannot be passed directly to the radix-5 decoder, because the place
/// values are different. The decoder reads the integer as radix-5 digits, so passing decimal `6`
/// directly would decode the radix-5 digit sequence `[1, 1]`, not the radix-4 digit sequence
/// `[1, 2]`.
///
/// The transformation therefore preserves the digit sequence, not the integer value. If a radix-4
/// motif index has digits:
///
/// `d0, d1, ..., d(k-1)`
///
/// where each digit is in `0..=3`, the helper builds the radix-5 code with the same digits:
///
/// `d0 * 5^(k-1) + d1 * 5^(k-2) + ... + d(k-1) * 5^0`
///
/// The loop does that left to right:
///
/// - `radix4_place` starts at `4^(k-1)`, the most significant radix-4 place.
/// - `digit = remaining / radix4_place` extracts the current radix-4 digit.
/// - `remaining %= radix4_place` removes that digit from the radix-4 value.
/// - `radix5_code = radix5_code * 5 + digit` appends the same digit to the radix-5 value.
/// - `radix4_place /= 4` moves to the next radix-4 place.
///
/// This is the same operation as decoding the radix-4 number into A/C/G/T digits and then encoding
/// those digits as a radix-5 k-mer code with no `N` digit.
///
/// For `k = 3`, radix-4 code `6` decomposes as:
///
/// `0 * 4^2 + 1 * 4^1 + 2 * 4^0`
///
/// Those digits are `[0, 1, 2]`, which mean `A`, `C`, and `G`. The matching radix-5 code is:
///
/// `0 * 5^2 + 1 * 5^1 + 2 * 5^0`
///
/// Passing that radix-5 code to `KmerSpec::decode_kmer()` produces `ACG`, while still using the
/// same decoder as the rest of the k-mer code. Iterating radix-4 codes from `0..4^k` therefore
/// produces the same order as iterating all radix-5 codes and filtering away motifs containing `N`,
/// but avoids decoding candidates that can never be part of the final motif axis.
///
/// Parameters
/// ----------
/// - `radix4_code`:
///   The 0-based index in the complete A/C/G/T-only motif set. It must be less than `4^k`.
/// - `k`:
///   Motif length in bases. `k = 0` maps to code `0`, representing the empty motif.
///
/// Returns
/// -------
/// - `radix5_code`:
///   The corresponding radix-5 code with no digit equal to `4`, ready for `KmerSpec::decode_kmer()`.
pub(crate) fn acgt_radix5_code_from_radix4(radix4_code: u64, k: usize) -> u64 {
    if k == 0 {
        return 0;
    }

    let mut remaining = radix4_code;
    let mut radix4_place = 4u64.pow((k - 1) as u32);
    let mut radix5_code = 0u64;

    for _ in 0..k {
        let digit = remaining / radix4_place;
        remaining %= radix4_place;
        radix5_code = radix5_code * 5 + digit;
        if radix4_place > 1 {
            radix4_place /= 4;
        }
    }

    radix5_code
}

/// Convert reduced selected k-mer counts into compact target-indexed output bins.
///
/// This is shared by commands whose motifs-file parser assigns original target indices before
/// counting starts. The caller supplies the sparse-weight validation and dense-output guard because
/// count semantics and dense-size limits belong to each command.
///
/// Unlike full-space k-mer postprocessing, this function does not decode motif strings from compact
/// k-mer keys. The motifs file has already defined the public target labels and assigned each
/// target an original index. This helper preserves that parser-assigned order for the public axis,
/// then remaps sparse row counts onto compact zero-based output columns.
///
/// Observed-only output keeps only targets that received at least one retained count. When
/// `include_all_targets` is true, every motifs-file target is retained and the caller-provided
/// dense-size guard runs before the dense writer is allowed to materialize the matrix.
///
/// Parameters
/// ----------
/// - `counts_by_window`:
///   Reduced sparse counts keyed by global output row and original motifs-file target index
/// - `total_windows`:
///   Number of rows in the final output matrix
/// - `labels`:
///   Motifs-file target labels in parser-assigned order
/// - `include_all_targets`:
///   Whether to keep unobserved motifs-file targets in the final output axis
/// - `should_store_weight`:
///   Command-specific sparse-weight validation and zeroish filter
/// - `ensure_dense_output_size`:
///   Command-specific dense-output size guard used when all targets are retained
///
/// Returns
/// -------
/// - `Result<(Vec<FxHashMap<u32, f64>>, Vec<String>)>`:
///   Compact output bins keyed by final column index, plus final column labels in order
pub(crate) fn postprocess_selected_motif_counts<RowCounts>(
    counts_by_window: impl IntoIterator<Item = (u64, RowCounts)>,
    total_windows: usize,
    labels: &[String],
    include_all_targets: bool,
    should_store_weight: fn(f64) -> Result<bool>,
    ensure_dense_output_size: fn(usize, usize) -> Result<()>,
) -> Result<(Vec<FxHashMap<u32, f64>>, Vec<String>)>
where
    RowCounts: IntoIterator<Item = (u32, f64)>,
{
    // Keep original target indices while deciding which motifs-file targets appear in output
    let mut raw_bins = vec![FxHashMap::default(); total_windows];
    let mut observed_targets = vec![false; labels.len()];

    for (original_idx, counts) in counts_by_window {
        // Validate the global row index before touching the output vector
        let idx: usize = original_idx
            .try_into()
            .context("selected motif window index does not fit in usize")?;
        ensure!(
            idx < raw_bins.len(),
            "selected motif window index {} is out of bounds for {} output windows",
            idx,
            raw_bins.len()
        );

        for (target_idx, value) in counts {
            // Tile merge already applies this check, but keeping it here protects the writer
            // boundary if another caller is added later
            if !should_store_weight(value)? {
                continue;
            }
            let target_position = target_idx as usize;
            ensure!(
                target_position < labels.len(),
                "selected motif target index {} is out of bounds for {} targets",
                target_position,
                labels.len()
            );
            observed_targets[target_position] = true;
            // Multiple tiles may contribute to the same final row and target
            *raw_bins[idx].entry(target_idx).or_insert(0.0) += value;
        }
    }

    // Build the public output axis in parser-assigned label order
    let mut output_index_by_target = vec![None; labels.len()];
    let mut motif_order = Vec::new();
    for (target_position, label) in labels.iter().enumerate() {
        if include_all_targets || observed_targets[target_position] {
            let output_index = u32::try_from(motif_order.len())
                .context("selected motif output index does not fit in u32")?;
            output_index_by_target[target_position] = Some(output_index);
            motif_order.push(label.clone());
        }
    }

    if include_all_targets {
        // All-target output may be dense, so fail before returning a matrix that cannot be written
        ensure_dense_output_size(total_windows, motif_order.len())?;
    }

    // Remap original motifs-file target indices to compact output column indices
    let mut indexed_bins = vec![FxHashMap::default(); total_windows];
    for (row_idx, raw_bin) in raw_bins.into_iter().enumerate() {
        for (target_idx, value) in raw_bin {
            let target_position = target_idx as usize;
            let output_index = output_index_by_target
                .get(target_position)
                .copied()
                .flatten()
                .with_context(|| {
                    format!(
                        "selected motif target index {} is missing from the output axis",
                        target_position
                    )
                })?;
            indexed_bins[row_idx].insert(output_index, value);
        }
    }

    Ok((indexed_bins, motif_order))
}

/// Split an aggregated `counts` map into per-k buckets.
///
/// * The `kmer_specs` dict tells us which k-values are valid and how to decode.
/// * Motifs that contain 'n' are discarded.
///
/// Returns one map for reference windows (“matches”) and one for mismatches.
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn split_and_decode_counts(
    counts: &FxHashMap<Kmer, f64>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
) -> DecodedCounts {
    let mut count_bins: FxHashMap<u8, FxHashMap<String, f64>> = FxHashMap::default();

    for (&kmer, &cnt) in counts {
        // User-readable motif, e.g. "ACG"
        let motif = kmer.to_string(kmer_specs);

        // Drop N's
        if motif.contains('N') {
            continue;
        }

        let bucket = count_bins.entry(kmer.k).or_default();
        *bucket.entry(motif).or_insert(0.0) += cnt;
    }

    DecodedCounts { counts: count_bins }
}

#[cfg(test)]
mod tests {
    include!("process_counts_tests.rs");
}
