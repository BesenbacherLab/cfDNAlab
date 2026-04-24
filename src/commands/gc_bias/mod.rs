pub mod binning;
#[cfg(feature = "cmd_gc_bias")]
pub mod config;
pub mod correct;
pub mod counting;
#[cfg(feature = "cmd_gc_bias")]
pub mod cross_tile_parts;
#[cfg(feature = "cmd_gc_bias")]
pub mod gc_bias;
pub mod interpolation;
pub mod load_reference_bias;
pub mod outliers;
pub mod package;
#[cfg(all(feature = "cmd_gc_bias", feature = "plotters"))]
pub mod plotting;
pub mod smoothing;
pub mod support_masking;
pub mod windows;

// Constants
pub const CORRECTION_CLAMP_RANGE: (f64, f64) = (0.1, 10.0);
pub const GC_CORRECTION_SCHEMA_VERSION: u32 = 2;
