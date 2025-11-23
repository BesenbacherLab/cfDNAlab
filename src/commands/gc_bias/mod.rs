pub mod binning;
pub mod config;
pub mod correct;
pub mod counting;
pub mod gc_bias;
pub mod interpolation;
pub mod load_reference_bias;
pub mod outliers;
pub mod package;
pub mod smoothing;
pub mod support_masking;

// Constants
pub const CORRECTION_CLAMP_RANGE: (f64, f64) = (0.1, 10.0);
pub const GC_CORRECTION_SCHEMA_VERSION: u32 = 1;
