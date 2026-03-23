use std::fmt::{Display, Formatter};

/// Error type for the public Rust API.
///
/// Use this to inspect why a public helper failed, for example when interval
/// bounds are invalid or an overlap fraction falls outside the supported range.
#[derive(Debug)]
pub enum Error {
    /// The interval bounds did not define a non-empty half-open interval.
    InvalidIntervalBounds { start: String, end: String },
    /// Converting an interval to another numeric coordinate type failed.
    InvalidIntervalConversion {
        start: String,
        end: String,
        target_type: &'static str,
    },
    /// Shifting an interval by the requested offset would overflow or underflow.
    InvalidIntervalOffset {
        start: String,
        end: String,
        offset: String,
    },
    /// The span bounds were inverted (`end < start`).
    InvalidSpanBounds { start: String, end: String },
    /// A fixed-bin helper received a bin size of zero.
    InvalidBinSize { bin_size: u64 },
    /// A fixed-size window index started at or past chromosome end.
    InvalidFixedWindowIndex {
        idx: u64,
        start: u64,
        chrom_len: u64,
    },
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
            Error::InvalidIntervalConversion {
                start,
                end,
                target_type,
            } => {
                write!(
                    formatter,
                    "converting interval [{start}, {end}) to {target_type} failed"
                )
            }
            Error::InvalidIntervalOffset { start, end, offset } => {
                write!(
                    formatter,
                    "offsetting interval [{start}, {end}) by {offset} would go out of bounds"
                )
            }
            Error::InvalidSpanBounds { start, end } => {
                write!(
                    formatter,
                    "span end ({end}) must be greater than or equal to start ({start})"
                )
            }
            Error::InvalidBinSize { bin_size } => {
                write!(formatter, "bin_size must be greater than 0, got {bin_size}")
            }
            Error::InvalidFixedWindowIndex {
                idx,
                start,
                chrom_len,
            } => {
                write!(
                    formatter,
                    "fixed-size window index {idx} starts at {start} beyond chromosome length {chrom_len}"
                )
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
