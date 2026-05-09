/// Minimum ACGT bases required when estimating GC fraction for sample reads.
pub const MIN_ACGT_BASES_FOR_GC_FRACTION: u32 = 10;

/// Maximum supported fragment length.
pub const MAX_SUPPORTED_FRAGMENT_LENGTH: u32 = 100_000;

/// Default maximum soft-clipped bases allowed at each fragment end.
pub const DEFAULT_MAX_SOFT_CLIPS: u16 = 256;

/// Maximum accepted value for per-end soft-clip expansion limits.
pub const MAX_MAX_SOFT_CLIPS: u16 = 256;

/// Version of the GC correction package schema.
pub const GC_CORRECTION_SCHEMA_VERSION: u32 = 2;

/// BAM AUX tag used for cfDNAlab-written GC correction weights.
pub const GC_WEIGHT_AUX_TAG: &[u8; 2] = b"GC";

/// BAM AUX tag used for cfDNAlab-written coverage-scaling weights.
pub const COVERAGE_WEIGHT_AUX_TAG: &[u8; 2] = b"cw";

/// BAM AUX tag used for cfDNAlab-written fragment-count-scaling weights.
pub const FRAGMENT_COUNT_WEIGHT_AUX_TAG: &[u8; 2] = b"nw";

/// BAM AUX tag used for cfDNAlab-written fragment lengths.
pub const FRAGMENT_LENGTH_AUX_TAG: &[u8; 2] = b"fl";
