use anyhow::{Result, ensure};
use std::fmt::Display;

/// A checked half-open interval `[start, end)`.
///
/// Use this for the geometric part of domain structs that carry genomic spans.
/// This type only represents non-empty intervals, so construction requires
/// `end > start`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Interval<T> {
    start: T,
    end: T,
}

impl<T> Interval<T>
where
    T: Copy + PartialOrd + Display,
{
    /// Create a checked half-open interval `[start, end)`.
    ///
    /// The interval must be non-empty, so `end` must be strictly greater than
    /// `start`.
    ///
    /// Use this when you want one place to enforce the half-open interval
    /// invariant instead of repeating start/end checks across callers.
    pub fn new(start: T, end: T) -> Result<Self> {
        ensure!(
            end > start,
            "interval end ({end}) must be greater than start ({start})"
        );
        Ok(Self { start, end })
    }

    /// Return the inclusive start coordinate.
    #[inline]
    pub fn start(&self) -> T {
        self.start
    }

    /// Return the exclusive end coordinate.
    #[inline]
    pub fn end(&self) -> T {
        self.end
    }

    /// Return the interval bounds as `(start, end)`.
    #[inline]
    pub fn into_inner(self) -> (T, T) {
        (self.start, self.end)
    }
}
