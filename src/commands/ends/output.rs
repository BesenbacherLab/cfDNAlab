//! Output-universe helpers for the `ends` command.
//!
//! This module owns the dense-vs-sparse motif ordering logic and the guards
//! around dense output size. Keeping that here avoids mixing output shape
//! policy into the main tile-processing runner.

use crate::{
    commands::ends::counting::format_end_motif_label,
    shared::{
        base::make_canonical,
        kmers::{kmer_codec::KmerSpec, process_counts::all_motifs as all_half_kmer_motifs},
    },
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;

const MAX_DENSE_END_MOTIF_OUTPUT_BYTES: u64 = 5 * 1024 * 1024 * 1024;

/// Collect the observed motif labels in deterministic sorted order.
///
/// This is the sparse-output path for `ends`. It keeps only motifs that actually
/// appeared in one or more windows and sorts them once at the end so downstream
/// output is stable across runs.
///
/// Parameters
/// ----------
/// - `bins`:
///   Final decoded per-window motif maps
///
/// Returns
/// -------
/// - `Vec<String>`:
///   Sorted observed motif labels
pub fn collect_end_motif_order(bins: &[FxHashMap<String, f64>]) -> Vec<String> {
    let mut motifs = std::collections::BTreeSet::new();
    for bin in bins {
        for motif in bin.keys() {
            motifs.insert(motif.clone());
        }
    }
    motifs.into_iter().collect()
}

/// Build the full dense motif universe for `--all-motifs`.
///
/// The returned labels follow the public `<outside>_<inside>` convention after
/// optional complement collapsing. This is only used for dense outputs where
/// every possible motif column must exist, even when a given sample has zero
/// counts for some motifs.
///
/// Parameters
/// ----------
/// - `inside_spec`:
///   Codec spec for the inside half, or `None` when `k_inside = 0`.
/// - `outside_spec`:
///   Codec spec for the outside half, or `None` when `k_outside = 0`.
/// - `collapse_complement`:
///   Whether complement-equivalent full motifs should collapse to one label.
///
/// Returns
/// -------
/// - `Result<Vec<String>>`:
///   Sorted dense motif universe for the current settings
pub fn build_all_end_motif_order(
    inside_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
    collapse_complement: bool,
) -> Result<Vec<String>> {
    let inside_motifs = all_half_motifs(inside_spec)?;
    let outside_motifs = all_half_motifs(outside_spec)?;
    let mut motifs = std::collections::BTreeSet::new();

    for outside in &outside_motifs {
        for inside in &inside_motifs {
            // Collapse operates on the biological motif sequence in `outside || inside` order.
            // The underscore is only a user-facing separator and is added after optional
            // same-orientation complement collapsing.
            let full_motif = format!("{outside}{inside}");
            let full_motif = if collapse_complement {
                make_canonical(full_motif, false, false)
            } else {
                full_motif
            };
            motifs.insert(format_end_motif_label(
                &full_motif,
                inside_spec,
                outside_spec,
            ));
        }
    }

    Ok(motifs.into_iter().collect())
}

/// Return all half-motif strings for one side of the end motif.
///
/// `None` means the corresponding `k` is zero, so the only valid half-motif is
/// the empty string.
///
/// Parameters
/// ----------
/// - `spec`:
///   Codec spec for one motif half, or `None` when that half is empty
///
/// Returns
/// -------
/// - `Result<Vec<String>>`:
///   All possible strings for that half
fn all_half_motifs(spec: Option<&KmerSpec>) -> Result<Vec<String>> {
    let Some(spec) = spec else {
        return Ok(vec![String::new()]);
    };

    let mut specs = FxHashMap::default();
    let k: u8 = spec
        .k
        .try_into()
        .context("k-mer size does not fit in u8 for motif enumeration")?;
    specs.insert(k, spec.clone());
    Ok(all_half_kmer_motifs(spec.k, &specs))
}

/// Guard dense output by the actual matrix size in bytes.
///
/// Dense end-motif output is convenient for small motif spaces because all
/// samples share the same columns. It becomes unreasonable once the matrix
/// itself would be too large, so this check fails early before allocation.
///
/// Parameters
/// ----------
/// - `n_windows`:
///   Number of output rows
/// - `n_motifs`:
///   Number of motif columns
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` when the dense matrix stays within the configured size budget
pub fn ensure_dense_end_motif_output_size(n_windows: usize, n_motifs: usize) -> Result<()> {
    let n_values = (n_windows as u64)
        .checked_mul(n_motifs as u64)
        .context("dense end-motif output shape overflows u64")?;
    let bytes = n_values
        .checked_mul(std::mem::size_of::<f64>() as u64)
        .context("dense end-motif output byte size overflows u64")?;

    ensure!(
        bytes <= MAX_DENSE_END_MOTIF_OUTPUT_BYTES,
        "Dense end-motif output would require {:.2} GiB for {} windows × {} motifs. \
         This output path is intentionally guarded. Reduce the motif space or window count for now.",
        bytes as f64 / (1024.0 * 1024.0 * 1024.0),
        n_windows,
        n_motifs
    );

    Ok(())
}

/// Guard full-universe enumeration for `--all-motifs`.
///
/// The dense path is defined by the final matrix size, not just by `k`, but
/// the full motif universe still needs its own early check so we do not try to
/// enumerate an obviously impossible motif list before the matrix is written.
///
/// Parameters
/// ----------
/// - `k_inside`:
///   Number of inside-fragment bases in the motif
/// - `k_outside`:
///   Number of outside-fragment bases in the motif
/// - `n_windows`:
///   Number of output rows
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` when dense `--all-motifs` enumeration stays within the size budget
pub fn ensure_all_motifs_enumeration_size(
    k_inside: usize,
    k_outside: usize,
    n_windows: usize,
) -> Result<()> {
    let total_k = k_inside
        .checked_add(k_outside)
        .context("combined motif length overflows usize")?;
    let motif_count_upper = 4_u64
        .checked_pow(total_k as u32)
        .context("all-motifs universe overflows u64")?;
    let n_motifs: usize = motif_count_upper
        .try_into()
        .context("all-motifs universe does not fit in usize")?;

    ensure_dense_end_motif_output_size(n_windows, n_motifs).with_context(|| {
        format!("refusing to enumerate all motifs for k_inside={k_inside}, k_outside={k_outside}")
    })?;

    Ok(())
}

#[cfg(test)]
mod tests {
    include!("output_tests.rs");
}
