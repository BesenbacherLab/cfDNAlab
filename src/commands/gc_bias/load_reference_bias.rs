use crate::{
    commands::cli_common::*,
    shared::bed::{Windows, load_scored_windows_from_bed},
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use ndarray::Array3;
use ndarray_npy::read_npy;
use std::path::Path;

pub struct ReferenceGcData {
    pub window_spec: WindowSpec,
    pub windows_map: Option<FxHashMap<String, Windows>>,
    pub window_indices_by_chr: Option<FxHashMap<String, Vec<u64>>>,
    pub counts: Array3<u64>,
    pub avg_window_size: Option<f64>,
}

pub fn load_reference_gc_data(
    ref_dir: &Path,
    chromosomes: Option<&[String]>,
    max_blacklisted_pct: u8,
) -> Result<ReferenceGcData> {
    let counts_path = ref_dir.join("all_ref_gc_counts.npy");
    let bins_path = ref_dir.join("ref_gc_bins.bed");

    let counts: Array3<u64> = read_npy(&counts_path)
        .with_context(|| format!("Reading reference GC counts from {:?}", counts_path))?;

    let window_spec = if counts.dim().0 == 1 && !bins_path.exists() {
        WindowSpec::Global
    } else {
        WindowSpec::Bed(bins_path.clone())
    };

    let (windows_map, avg_window_size) = if matches!(window_spec, WindowSpec::Bed(_)) {
        let mut windows_map = parse_reference_bins(&bins_path, chromosomes, max_blacklisted_pct)
            .with_context(|| format!("Reading reference GC window coordinates {:?}", bins_path))?;

        // Ensure we keep a Windows object per requested chromosome, even when empty
        if let Some(chroms) = chromosomes {
            for chr in chroms {
                windows_map
                    .entry(chr.clone())
                    .or_insert_with(|| Windows::from_sorted(Vec::new()));
            }
        }

        // Windows are filtered by blacklist percentage
        let total_windows: usize = windows_map.iter().map(|(_, ws)| ws.len()).sum();
        ensure!(
            total_windows > 0,
            "Reference GC BED does not provide any windows. Please supply a non-empty BED file"
        );

        let sum_of_spans = windows_map
            .iter()
            .map(|(_, ws)| ws.as_slice().iter().map(|w| w.1 - w.0).sum::<u64>())
            .sum::<u64>();
        let avg_window_span = sum_of_spans as f64 / total_windows as f64;
        ensure!(
            counts.shape()[0] <= total_windows,
            "Found more reference window coordinates than count distribution windows",
        );

        (Some(windows_map), Some(avg_window_span))
    } else {
        (None, None)
    };

    let window_indices = if let Some(ws_map) = &windows_map {
        let window_indices: FxHashMap<String, Vec<u64>> = ws_map
            .iter()
            .map(|(chr, ws)| {
                (
                    chr.to_owned(),
                    ws.windows.iter().map(|(_, _, idx)| *idx).collect(),
                )
            })
            .collect();
        Some(window_indices)
    } else {
        None
    };

    Ok(ReferenceGcData {
        window_spec: window_spec,
        windows_map: windows_map,
        window_indices_by_chr: window_indices,
        counts,
        avg_window_size: avg_window_size,
    })
}

pub fn parse_reference_bins(
    bins_path: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    max_blacklisted_pct: u8,
) -> Result<FxHashMap<String, Windows>> {
    let bins_path = bins_path.as_ref();

    // Filter function for removing windows with high blacklisting
    let threshold = max_blacklisted_pct as f64 / 100.0;
    let filter_windows_fn: &dyn Fn(&str, u64, u64, f64) -> bool =
        &move |_: &str, _: u64, _: u64, pct: f64| pct > threshold;

    let scored_windows_map: FxHashMap<String, Windows> =
        load_scored_windows_from_bed(bins_path, chromosomes, Some(filter_windows_fn))?
            .iter()
            .map(|(chr, ws)| (chr.to_owned(), ws.to_windows()))
            .collect();

    Ok(scored_windows_map)
}
