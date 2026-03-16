use std::fmt::{Display, Formatter};

/// Error type for the public Rust API.
///
/// Use this to inspect why a public helper failed, for example when interval
/// bounds are invalid or an overlap fraction falls outside the supported range.
#[derive(Debug)]
pub enum Error {
    /// The interval bounds did not define a non-empty half-open interval.
    InvalidIntervalBounds { start: String, end: String },
    /// A fixed-bin helper received a bin size of zero.
    InvalidBinSize { bin_size: u64 },
    /// An overlap fraction fell outside the inclusive range `[0.0, 1.0]`.
    OverlapFractionOutOfBounds { overlap_fraction: f32 },
    /// A tile fetch range did not fully cover the tile core.
    TileFetchDoesNotCoverCore,
}

impl Display for Error {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidIntervalBounds { start, end } => {
                write!(
                    formatter,
                    "interval end ({end}) must be greater than start ({start})"
                )
            }
            Error::InvalidBinSize { bin_size } => {
                write!(formatter, "bin_size must be greater than 0, got {bin_size}")
            }
            Error::OverlapFractionOutOfBounds { overlap_fraction } => {
                write!(
                    formatter,
                    "overlap_fraction was out of bounds (0.0-1.0): {overlap_fraction}"
                )
            }
            Error::TileFetchDoesNotCoverCore => {
                write!(
                    formatter,
                    "tile fetch interval must fully cover the tile core"
                )
            }
        }
    }
}

impl std::error::Error for Error {}

/// Result alias for the public Rust API.
pub type Result<T> = std::result::Result<T, Error>;
