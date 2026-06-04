use crate::commands::ref_gc_bias::zarr::read_reference_gc_package_zarr;
use crate::shared::reference::ContigFootprintEntry;
use anyhow::Result;
use ndarray::Array2;
use std::path::Path;

#[derive(Clone, Debug)]
pub(crate) struct ReferenceGCData {
    pub(crate) counts: Array2<f64>,
    pub(crate) unobservables_support_mask: Array2<bool>,
    pub(crate) outliers_support_mask: Array2<bool>,
    pub(crate) gc_percent_widths: Array2<u16>,
    pub(crate) metadata: ReferenceGCMetadata,
}

#[derive(Clone, Debug)]
pub(crate) struct ReferenceGCMetadata {
    pub(crate) min_fragment_length: usize,
    pub(crate) max_fragment_length: usize,
    pub(crate) end_offset: u8,
    pub(crate) chromosomes: Vec<String>,
    pub(crate) reference_contig_footprint: Vec<ContigFootprintEntry>,
    pub(crate) skip_interpolation: bool,
    pub(crate) smoothing_sigma: f64,
    pub(crate) smoothing_radius: u8,
    pub(crate) skip_smoothing: bool,
}

pub(crate) fn load_reference_gc_data(ref_file: &Path) -> Result<ReferenceGCData> {
    read_reference_gc_package_zarr(ref_file)
}
