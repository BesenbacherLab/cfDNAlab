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
