use std::{fs::File, path::Path};

use anyhow::{Context, Result, ensure};
use ndarray::{Array1, Array2, ArrayView1};
use ndarray_npy::{NpzReader, NpzWriter, ReadNpyExt};

use crate::{
    commands::{cli_common::WindowSpec, lengths::counting::LengthCounts},
    shared::tiled_run::{
        Tile, TileWindowSpan, clamp_fetch_to_window_span, parse_tile_index, tile_window_min_max,
    },
};

/// Write per-tile partial length counts as an NPZ archive.
///
/// The archive stores two arrays:
/// - `window_idx_chr` (u64): Zero-based window index within the current chromosome.
/// - `counts` (f64): Matrix with one row per window and one column per fragment length.
pub fn write_partials_npz(
    temp_dir: &Path,
    prefix: &str,
    chr: &str,
    tile_idx: u32,
    window_idxs_chr: &[u64],
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
            "All count rows must have identical length"
        );
    }

    let path = temp_dir.join(format!("{prefix}.{chr}.{tile_idx}.npz"));
    let file = File::create(&path)?;
    let mut npz = NpzWriter::new(file);

    let idx_arr = Array1::from(window_idxs_chr.to_vec());
    let mut counts_arr = Array2::zeros((counts.len(), counts_len));
    for (row_idx, lc) in counts.iter().enumerate() {
        let row = ArrayView1::from(lc.counts.as_slice());
        counts_arr.row_mut(row_idx).assign(&row);
    }

    npz.add_array("window_idx_chr", &idx_arr)?;
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
    temp_dir: &Path,
    partials_prefix: &str,
    cross_prefix: &str,
    n_windows: usize,
    template: &LengthCounts,
) -> Result<Vec<LengthCounts>> {
    // Expected contributions per window, incremented by crossing files and set to 1 later for contained windows
    let mut expected: Vec<u32> = vec![0; n_windows];
    // Actual contributions observed while merging partials
    let mut contributions: Vec<u32> = vec![0; n_windows];
    // Accumulator for summed counts per window
    let mut counts_by_idx: Vec<LengthCounts> = std::iter::repeat_with(|| template.zeroed_like())
        .take(n_windows)
        .collect();

    // First accumulate expected contributions from crossing windows
    for entry in
        std::fs::read_dir(temp_dir).with_context(|| format!("Listing {}", temp_dir.display()))?
    {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let fname = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
        if !fname.starts_with(cross_prefix) || !fname.contains(&format!(".{chr}.")) {
            continue;
        }
        // Skip files that do not carry a numeric tile index suffix
        if parse_tile_index(fname).is_none() {
            continue;
        }

        let file =
            File::open(&path).with_context(|| format!("Opening cross file {}", path.display()))?;
        let arr: Array1<u64> = Array1::read_npy(file)
            .with_context(|| format!("Reading cross file {}", path.display()))?;
        for idx in arr.iter() {
            let i = *idx as usize;
            ensure!(
                i < expected.len(),
                "Cross index {} out of bounds for chromosome {}",
                idx,
                chr
            );
            expected[i] = expected[i].saturating_add(1);
        }
    }

    // Windows never listed in crossing files are contained within a single tile, so they expect one contribution
    for exp in expected.iter_mut() {
        if *exp == 0 {
            *exp = 1;
        }
    }

    // Then sum partial counts
    for entry in
        std::fs::read_dir(temp_dir).with_context(|| format!("Listing {}", temp_dir.display()))?
    {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let fname = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
        if !fname.starts_with(partials_prefix) || !fname.contains(&format!(".{chr}.")) {
            continue;
        }
        // Skip files that do not carry a numeric tile index suffix
        if parse_tile_index(fname).is_none() {
            continue;
        }

        let file =
            File::open(&path).with_context(|| format!("Opening partials {}", path.display()))?;
        let mut npz =
            NpzReader::new(file).with_context(|| format!("Reading partials {}", path.display()))?;
        let idxs: Array1<u64> = npz
            .by_name("window_idx_chr")
            .with_context(|| format!("Reading window_idx_chr in {}", path.display()))?;
        let counts: Array2<f64> = npz
            .by_name("counts")
            .with_context(|| format!("Reading counts in {}", path.display()))?;
        ensure!(
            counts.nrows() == idxs.len(),
            "counts rows did not match idx length in {}",
            path.display()
        );
        ensure!(
            counts.ncols() == template.counts.len(),
            "counts width mismatch in {}",
            path.display()
        );

        for (row_idx, idx) in idxs.iter().enumerate() {
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
        }
    }

    // Validate contributions
    for (i, (have, want)) in contributions.iter().zip(expected.iter()).enumerate() {
        ensure!(
            *have == *want,
            "Window {} on {} had {} contributions but expected {}",
            i,
            chr,
            have,
            want
        );
    }

    ensure!(
        contributions.iter().all(|c| *c > 0),
        "Some windows received zero contributions on {}",
        chr
    );

    Ok(counts_by_idx)
}

/// Determine the fetch span for a tile based on the active window strategy.
///
/// Global mode: the full tile fetch range is used.
///
/// Fixed-size window mode: the span is defined as the first and last bin that
/// touches the tile core, clamped to the chromosome length.
///
/// BED mode: uses the precomputed window bounds for the tile to avoid fetching
/// unrelated regions and returns `None` when the tile does not intersect any
/// BED windows.
///
/// Parameters
/// ----------
///
/// - `tile`: Tile describing the chromosome, core span, and fetch span.
///
/// - `tile_window_span`: Cached min and max window bounds for the tile (BED mode).
///
/// - `windows_chr`: Chromosome BED windows as `(start, end, idx)` tuples (BED mode).
///
/// - `window_opt`: Window specification selecting global, fixed-size, or BED mode.
///
/// - `chrom_len`: Chromosome length used to clamp fetch coordinates.
///
/// Returns
/// - `(start, end)` fetch coordinates as `i64`, or `None` when no windows apply.
pub fn fetch_span_for_tile(
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_chr: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    chrom_len: u64,
) -> Option<(i64, i64)> {
    match window_opt {
        WindowSpec::Global => Some((
            tile.fetch_start as i64,
            (tile.fetch_end.min(chrom_len as u32)) as i64,
        )),
        WindowSpec::Size(window_bp) => {
            let core_start = tile.core_start as u64;
            let core_end = (tile.core_end as u64).min(chrom_len);
            if core_start >= chrom_len || core_end == 0 {
                return None;
            }
            let window_idx_start = core_start / window_bp;
            let window_idx_end = (core_end.saturating_sub(1)) / window_bp;
            let window_start = window_idx_start * window_bp;
            let window_end = ((window_idx_end + 1) * window_bp).min(chrom_len);
            clamp_fetch_to_window_span(tile, chrom_len, window_start, window_end)
        }
        WindowSpec::Bed(_) => {
            let wchr = windows_chr?;
            let (min_ws, max_we) = tile_window_min_max(wchr, tile, tile_window_span)?;
            clamp_fetch_to_window_span(tile, chrom_len, min_ws, max_we)
        }
    }
}
