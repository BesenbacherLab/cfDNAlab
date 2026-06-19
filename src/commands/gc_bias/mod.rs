pub(crate) mod binning;
pub(crate) mod config;
pub(crate) mod correct;
pub(crate) mod counting;
pub(crate) mod cross_tile_parts;
pub(crate) mod gc_bias;
pub(crate) mod interpolation;
pub(crate) mod load_reference_bias;
pub(crate) mod outliers;
pub(crate) mod package;
#[cfg(feature = "plotters")]
pub(crate) mod plotting;
pub(crate) mod smoothing;
pub(crate) mod support_masking;
pub(crate) mod windows;

// Constants
pub(crate) const CORRECTION_CLAMP_RANGE: (f64, f64) = (0.1, 10.0);
