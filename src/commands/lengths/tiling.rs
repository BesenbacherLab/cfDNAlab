use std::{
    fs::File,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, ensure};
use ndarray::{Array1, Array2, ArrayView1};
use ndarray_npy::{NpzReader, NpzWriter, ReadNpyExt};

use crate::commands::lengths::counting::LengthCounts;

/// Write per-tile partial length counts as an NPZ archive.
///
/// The archive stores two arrays:
/// - `window_idx_chr` (u64): Zero-based window index within the current chromosome.
/// - `counts` (f64): Matrix with one row per window and one column per
///   fragment length.
pub fn write_partials_npz(
    temp_dir: &Path,
    prefix: &str,
    chr: &str,
    tile_idx: u32,
    window_idxs_chr: &[u64],
    contained_flags: &[bool],
    counts: &[LengthCounts],
) -> Result<Option<std::path::PathBuf>> {
    if window_idxs_chr.is_empty() {
        return Ok(None);
    }

    let counts_len = counts
        .first()
        .map(|c| c.counts.len())
        .context("counts array empty")?;
    for c in counts {
        ensure!(
            c.counts.len() == counts_len,
            "All count vectors (rows) must have identical length"
        );
    }

    ensure!(
        contained_flags.len() == window_idxs_chr.len(),
        "contained flags length mismatch for tile {} {}",
        chr,
        tile_idx
    );
    ensure!(
        counts.len() == window_idxs_chr.len(),
        "counts length mismatch for tile {} {}",
        chr,
        tile_idx
    );

    let path = temp_dir.join(format!("{prefix}.{chr}.{tile_idx}.npz"));
    let file = File::create(&path)?;
    let mut npz = NpzWriter::new(file);

    let idx_arr = Array1::from(window_idxs_chr.to_vec());
    let contained_arr: Array1<u8> = contained_flags
        .iter()
        .map(|&b| if b { 1u8 } else { 0u8 })
        .collect();
    let mut counts_arr = Array2::zeros((counts.len(), counts_len));
    for (row_idx, lc) in counts.iter().enumerate() {
        let row = ArrayView1::from(lc.counts.as_slice());
        counts_arr.row_mut(row_idx).assign(&row);
    }

    npz.add_array("window_idx_chr", &idx_arr)?;
    npz.add_array("contained", &contained_arr)?;
    npz.add_array("counts", &counts_arr)?;
    npz.finish()?;
    Ok(Some(path))
}

/// Write the list of crossing window indices as an NPY array.
pub fn write_cross_npy(
    temp_dir: &Path,
    prefix: &str,
    chr: &str,
    tile_idx: u32,
    crossing_window_idxs_chr: &[u64],
) -> Result<Option<std::path::PathBuf>> {
    if crossing_window_idxs_chr.is_empty() {
        return Ok(None);
    }
    let path = temp_dir.join(format!("{prefix}.{chr}.{tile_idx}.npy"));
    let arr = Array1::from(crossing_window_idxs_chr.to_vec());
    ndarray_npy::write_npy(&path, &arr)?;
    Ok(Some(path))
}

/// Reduce partial NPZ files and crossing NPY files for one chromosome into full-length counts.
///
/// `n_windows` must equal the total number of window indices for the chromosome (scan order).
pub fn reduce_partials_for_chr(
    chr: &str,
    partial_paths: &[PathBuf],
    cross_paths: &[PathBuf],
    n_windows: usize,
    template: &LengthCounts,
) -> Result<Vec<LengthCounts>> {
    // Expected contributions per window, tracked separately for crossing and contained tiles
    let mut cross_counts: Vec<u32> = vec![0; n_windows];
    let mut contained_counts: Vec<u32> = vec![0; n_windows];
    let mut contributions: Vec<u32> = vec![0; n_windows];
    // Accumulator for summed counts per window
    let mut counts_by_idx: Vec<LengthCounts> = std::iter::repeat_with(|| template.zeroed_like())
        .take(n_windows)
        .collect();

    // First accumulate contributions from crossing windows
    for path in cross_paths {
        let file =
            File::open(path).with_context(|| format!("Opening cross file {}", path.display()))?;
        let arr: Array1<u64> = Array1::read_npy(file)
            .with_context(|| format!("Reading cross file {}", path.display()))?;
        for idx in arr.iter() {
            let i = *idx as usize;
            ensure!(
                i < cross_counts.len(),
                "Cross index {} out of bounds for chromosome {}",
                idx,
                chr
            );
            cross_counts[i] = cross_counts[i].saturating_add(1);
        }
    }

    // Then sum partial counts
    for path in partial_paths {
        let file =
            File::open(path).with_context(|| format!("Opening partials {}", path.display()))?;
        let mut npz =
            NpzReader::new(file).with_context(|| format!("Reading partials {}", path.display()))?;
        let idxs: Array1<u64> = npz
            .by_name("window_idx_chr")
            .with_context(|| format!("Reading window_idx_chr in {}", path.display()))?;
        let contained: Array1<u8> = npz
            .by_name("contained")
            .with_context(|| format!("Reading contained in {}", path.display()))?;
        let counts: Array2<f64> = npz
            .by_name("counts")
            .with_context(|| format!("Reading counts in {}", path.display()))?;
        ensure!(
            counts.nrows() == idxs.len(),
            "counts rows did not match idx length in {}",
            path.display()
        );
        ensure!(
            contained.len() == idxs.len(),
            "contained length did not match idx length in {}",
            path.display()
        );
        ensure!(
            counts.ncols() == template.counts.len(),
            "counts width mismatch in {}",
            path.display()
        );

        for (row_idx, (idx, contained_flag)) in idxs.iter().zip(contained.iter()).enumerate() {
            let i = *idx as usize;
            ensure!(
                i < counts_by_idx.len(),
                "Partial index {} out of bounds for chromosome {}",
                idx,
                chr
            );
            let row_view = counts.row(row_idx);
            let row_slice = row_view.as_slice().context("counts row not contiguous")?;
            for (dst, val) in counts_by_idx[i].counts.iter_mut().zip(row_slice.iter()) {
                *dst += *val;
            }
            contributions[i] = contributions[i].saturating_add(1);
            if *contained_flag == 1 {
                contained_counts[i] = contained_counts[i].saturating_add(1);
            }
        }
    }

    // Validate contributions
    for (i, have) in contributions.iter().enumerate() {
        let expected = cross_counts[i].saturating_add(contained_counts[i]).max(1);
        ensure!(
            *have == expected,
            "Window {} on {} had {} contributions but expected {} (cross files counted {}, contained {})",
            i,
            chr,
            have,
            expected,
            cross_counts[i],
            contained_counts[i]
        );
    }

    ensure!(
        contributions.iter().all(|c| *c > 0),
        "Some windows received zero contributions on {}",
        chr
    );

    Ok(counts_by_idx)
}
