use crate::commands::gc_bias::GC_CORRECTION_SCHEMA_VERSION;
use anyhow::{Context, Result, ensure};
use ndarray::{Array1, Array2};
use ndarray_npy::NpzReader;
use std::{fs::File, path::Path};

#[derive(Clone, Debug)]
pub struct ReferenceGCData {
    pub counts: Array2<f64>,
    pub unobservables_support_mask: Array2<bool>,
    pub outliers_support_mask: Array2<bool>,
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

pub fn load_reference_gc_data(ref_file: &Path) -> Result<ReferenceGCData> {
    let (counts, unobservables_support_mask, outliers_support_mask, gc_percent_widths, metadata) =
        read_reference_gc_package(&ref_file)?;

    Ok(ReferenceGCData {
        counts,
        unobservables_support_mask,
        outliers_support_mask,
        gc_percent_widths,
        metadata,
    })
}

fn read_reference_gc_package(
    path: &Path,
) -> Result<(
    Array2<f64>,
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

    let counts: Array2<f64> = reader
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
    let version_arr: Array1<u32> = reader
        .by_name("version")
        .context("missing version in reference GC package")?;
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
        version_arr.len() == 1,
        "version should be length 1. Found len={}",
        version_arr.len()
    );
    ensure!(
        end_offset_arr.len() == 1,
        "end_offset should be length 1. Found len={}",
        end_offset_arr.len()
    );
    ensure!(
        skip_interpolation_arr.len() == 1,
        "skip_interpolation should be length 1. Found len={}",
        skip_interpolation_arr.len()
    );
    ensure!(
        smoothing_radius_arr.len() == 1,
        "smoothing_radius should be length 1. Found len={}",
        smoothing_radius_arr.len()
    );
    ensure!(
        smoothing_sigma_arr.len() == 1,
        "smoothing_sigma should be length 1. Found len={}",
        smoothing_sigma_arr.len()
    );
    ensure!(
        skip_smoothing_arr.len() == 1,
        "skip_smoothing should be length 1. Found len={}",
        skip_smoothing_arr.len()
    );
    ensure!(
        version_arr[0] == GC_CORRECTION_SCHEMA_VERSION,
        "Reference GC package schema version mismatch: file={}, expected={}; \
        Incompatible with this version of cfDNAlab.",
        version_arr[0],
        GC_CORRECTION_SCHEMA_VERSION
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
        unobservables_support_mask.dim() == counts.dim(),
        "Reference counts ({:?}) and support masks ({:?}) had incompatible shapes",
        counts.dim(),
        unobservables_support_mask.dim()
    );
    ensure!(
        gc_percent_widths.dim() == counts.dim(),
        "GC percent widths shape {:?} must match per-window counts shape {:?}",
        gc_percent_widths.dim(),
        counts.dim(),
    );

    Ok((
        counts,
        unobservables_support_mask,
        outliers_support_mask,
        gc_percent_widths,
        metadata,
    ))
}
