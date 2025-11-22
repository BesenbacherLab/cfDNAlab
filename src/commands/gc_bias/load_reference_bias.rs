use crate::{
    commands::cli_common::*,
    shared::bed::{Windows, load_scored_windows_from_bed},
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use ndarray::{Array2, Array3};
use ndarray_npy::read_npy;
use std::path::Path;

pub struct ReferenceGCData {
    pub window_spec: WindowSpec,
    pub windows_map: Option<FxHashMap<String, Windows>>,
    pub window_indices_by_chr: Option<FxHashMap<String, Vec<u64>>>,
    pub counts: Array3<f64>,
    pub unobservables_support_mask: Array2<bool>,
    pub outliers_support_mask: Array2<bool>,
    pub avg_window_size: Option<f64>,
}

pub fn load_reference_gc_data(
    ref_dir: &Path,
    chromosomes: Option<&[String]>,
    max_blacklisted_pct: u8,
) -> Result<ReferenceGCData> {
    let counts_path = ref_dir.join("ref_gc_counts.npy");
    let unobservables_support_mask_path = ref_dir.join("ref_support_mask.unobservables.npy");
    let outliers_support_mask_path = ref_dir.join("ref_support_mask.outliers.npy");
    let bins_path = ref_dir.join("ref_gc_bins.bed");

    let counts: Array3<f64> = read_npy(&counts_path)
        .with_context(|| format!("Reading reference GC counts from {:?}", counts_path))?;

    let unobservables_support_mask: Array2<bool> = read_npy(&unobservables_support_mask_path)
        .with_context(|| {
            format!(
                "Reading reference support mask (unobservables) from {:?}",
                unobservables_support_mask_path
            )
        })?;

    let outliers_support_mask: Array2<bool> =
        read_npy(&outliers_support_mask_path).with_context(|| {
            format!(
                "Reading reference support mask (outliers) from {:?}",
                outliers_support_mask_path
            )
        })?;

    ensure!(
        unobservables_support_mask.dim() == outliers_support_mask.dim(),
        "The two support masks must have the same shape. Unobservables ({:?}) != outliers ({:?}).",
        unobservables_support_mask.dim(),
        outliers_support_mask.dim(),
    );

    let num_count_windows = counts.dim().0;

    ensure!(
        unobservables_support_mask.dim().0 == counts.dim().1
            && unobservables_support_mask.dim().1 == counts.dim().2,
        "Reference counts ({:?}) and support masks ({:?}) had incompatible shapes",
        counts.dim(),
        unobservables_support_mask.dim()
    );

    let window_spec = if num_count_windows == 1 && !bins_path.exists() {
        WindowSpec::Global
    } else {
        WindowSpec::Bed(bins_path.clone())
    };

    let (windows_map, avg_window_size) = if matches!(window_spec, WindowSpec::Bed(_)) {
        let mut windows_map = parse_reference_bins(
            &bins_path,
            chromosomes,
            max_blacklisted_pct,
            num_count_windows as u64,
        )
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

    Ok(ReferenceGCData {
        window_spec: window_spec,
        windows_map: windows_map,
        window_indices_by_chr: window_indices,
        counts,
        unobservables_support_mask,
        outliers_support_mask,
        avg_window_size: avg_window_size,
    })
}

pub fn parse_reference_bins(
    bins_path: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    max_blacklisted_pct: u8,
    exp_num_windows: u64,
) -> Result<FxHashMap<String, Windows>> {
    let bins_path = bins_path.as_ref();

    // Filter function for removing windows with high blacklisting
    let threshold = max_blacklisted_pct as f64 / 100.0;
    let filter_windows_fn: &dyn Fn(&str, u64, u64, f64) -> bool =
        &move |_: &str, _: u64, _: u64, pct: f64| pct <= threshold;

    let scored_windows_map: FxHashMap<String, Windows> = load_scored_windows_from_bed(
        bins_path,
        chromosomes,
        Some(filter_windows_fn),
        Some(exp_num_windows),
    )?
    .iter()
    .map(|(chr, ws)| (chr.to_owned(), ws.to_windows()))
    .collect();

    Ok(scored_windows_map)
}
