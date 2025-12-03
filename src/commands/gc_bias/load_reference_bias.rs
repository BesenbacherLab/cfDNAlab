use crate::{
    commands::cli_common::*,
    shared::bed::{Windows, load_scored_windows_from_bed},
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use ndarray::{Array1, Array2, Array3};
use ndarray_npy::NpzReader;
use std::{fs::File, path::Path};

#[derive(Clone, Debug)]
pub struct ReferenceGCData {
    pub window_spec: WindowSpec,
    pub windows_map: Option<FxHashMap<String, Windows>>,
    pub window_indices_by_chr: Option<FxHashMap<String, Vec<u64>>>,
    pub counts: Array3<f64>,
    pub unobservables_support_mask: Array2<bool>,
    pub outliers_support_mask: Array2<bool>,
    pub avg_window_size: Option<f64>,
    pub gc_percent_widths: Array2<u16>,
    pub metadata: ReferenceGCMetadata,
}

#[derive(Clone, Debug)]
pub struct ReferenceGCMetadata {
    pub min_fragment_length: usize,
    pub max_fragment_length: usize,
    pub end_offset: u8,
    pub skip_interpolation: bool,
    pub smoothing_sigma: f64,
    pub smoothing_radius: u8,
    pub skip_smoothing: bool,
}

pub fn load_reference_gc_data(
    ref_dir: &Path,
    chromosomes: Option<&[String]>,
    max_blacklisted_pct: u8,
) -> Result<ReferenceGCData> {
    let package_path = ref_dir.join("ref_gc_package.npz");
    let bins_path = ref_dir.join("ref_gc_bins.bed");
    let (counts, unobservables_support_mask, outliers_support_mask, gc_percent_widths, metadata) =
        read_reference_gc_package(&package_path)?;

    let num_count_windows = counts.dim().0;

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
        gc_percent_widths,
        metadata,
    })
}

fn read_reference_gc_package(
    path: &Path,
) -> Result<(
    Array3<f64>,
    Array2<bool>,
    Array2<bool>,
    Array2<u16>,
    ReferenceGCMetadata,
)> {
    let file = File::open(path).with_context(|| {
        format!(
            "Reading reference GC package from {:?}. Regenerate reference GC if missing",
            path
        )
    })?;
    let mut reader = NpzReader::new(file)?;

    let counts: Array3<f64> = reader
        .by_name("counts")
        .context("Missing counts in reference GC package")?;
    let unobservables_support_mask: Array2<bool> = reader
        .by_name("support_mask_unobservables")
        .context("Missing support_mask_unobservables in reference GC package")?;
    let outliers_support_mask: Array2<bool> = reader
        .by_name("support_mask_outliers")
        .context("Missing support_mask_outliers in reference GC package")?;
    let gc_percent_widths: Array2<u16> = reader
        .by_name("gc_percent_widths")
        .context("Missing gc_percent_widths in reference GC package")?;
    let lengths: Array1<u32> = reader
        .by_name("length_range")
        .context("missing length_range in reference GC package")?;
    ensure!(
        lengths.len() == 2,
        "length_range should contain [min, max] (len=2). Found len={}",
        lengths.len()
    );
    let end_offset_arr: Array1<u32> = reader
        .by_name("end_offset")
        .context("missing end_offset in reference GC package")?;
    let skip_interpolation_arr: Array1<bool> = reader
        .by_name("skip_interpolation")
        .context("missing skip_interpolation in reference GC package")?;
    let smoothing_radius_arr: Array1<u32> = reader
        .by_name("smoothing_radius")
        .context("missing smoothing_radius in reference GC package")?;
    let smoothing_sigma_arr: Array1<f64> = reader
        .by_name("smoothing_sigma")
        .context("missing smoothing_sigma in reference GC package")?;
    let skip_smoothing_arr: Array1<bool> = reader
        .by_name("skip_smoothing")
        .context("missing skip_smoothing in reference GC package")?;
    ensure!(
        end_offset_arr.len() == 1,
        "end_offset should be length 1. Found len={}",
        end_offset_arr.len()
    );
    let metadata = ReferenceGCMetadata {
        min_fragment_length: lengths[0] as usize,
        max_fragment_length: lengths[1] as usize,
        end_offset: end_offset_arr[0] as u8,
        skip_interpolation: skip_interpolation_arr[0] as bool,
        smoothing_radius: smoothing_radius_arr[0] as u8,
        smoothing_sigma: smoothing_sigma_arr[0] as f64,
        skip_smoothing: skip_smoothing_arr[0] as bool,
    };

    ensure!(
        unobservables_support_mask.dim() == outliers_support_mask.dim(),
        "The two support masks must have the same shape. Unobservables ({:?}) != outliers ({:?}).",
        unobservables_support_mask.dim(),
        outliers_support_mask.dim(),
    );
    ensure!(
        unobservables_support_mask.dim().0 == counts.dim().1
            && unobservables_support_mask.dim().1 == counts.dim().2,
        "Reference counts ({:?}) and support masks ({:?}) had incompatible shapes",
        counts.dim(),
        unobservables_support_mask.dim()
    );
    ensure!(
        gc_percent_widths.dim() == (counts.dim().1, counts.dim().2),
        "GC percent widths shape {:?} must match per-window counts shape {:?}",
        gc_percent_widths.dim(),
        (counts.dim().1, counts.dim().2)
    );

    Ok((
        counts,
        unobservables_support_mask,
        outliers_support_mask,
        gc_percent_widths,
        metadata,
    ))
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
