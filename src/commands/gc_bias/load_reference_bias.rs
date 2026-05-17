use crate::commands::ref_gc_bias::zarr::read_reference_gc_package_zarr;
use crate::shared::reference::ContigFootprintEntry;
use anyhow::Result;
use ndarray::Array2;
use std::path::Path;

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
    pub chromosomes: Vec<String>,
    pub reference_contig_footprint: Vec<ContigFootprintEntry>,
    pub skip_interpolation: bool,
    pub smoothing_sigma: f64,
    pub smoothing_radius: u8,
    pub skip_smoothing: bool,
}

pub fn load_reference_gc_data(ref_file: &Path) -> Result<ReferenceGCData> {
    read_reference_gc_package_zarr(ref_file)
}
