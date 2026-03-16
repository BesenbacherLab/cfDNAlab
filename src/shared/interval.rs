use crate::{Error, Result};
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
        if end <= start {
            return Err(Error::InvalidIntervalBounds {
                start: start.to_string(),
                end: end.to_string(),
            });
        }
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

impl<T> TryFrom<(T, T)> for Interval<T>
where
    T: Copy + PartialOrd + Display,
{
    type Error = Error;

    /// Convert a `(start, end)` tuple into a checked half-open interval.
    ///
    /// Use this when coordinates are already stored as tuples and you want to
    /// validate them without unpacking them manually. This is especially useful
    /// when collecting many tuples into `Vec<Interval<_>>`.
    ///
    /// Parameters
    /// ----------
    /// - `bounds`:
    ///   Interval bounds as `(start, end)`.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Checked non-empty interval.
    ///
    /// Example
    /// -------
    /// ```rust
    /// use cfdnalab::shared::interval::Interval;
    ///
    /// let intervals: cfdnalab::Result<Vec<_>> = vec![(5_u64, 6_u64), (10, 20)]
    ///     .into_iter()
    ///     .map(Interval::try_from)
    ///     .collect();
    ///
    /// assert_eq!(intervals?.len(), 2);
    /// # Ok::<(), cfdnalab::Error>(())
    /// ```
    fn try_from(bounds: (T, T)) -> Result<Self> {
        Self::new(bounds.0, bounds.1)
    }
}
