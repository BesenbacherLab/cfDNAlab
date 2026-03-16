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

/// A checked half-open interval together with an external index or identifier.
///
/// Use this when an interval needs stable caller metadata, such as the original
/// window index from a BED file. The interval part still follows the same
/// non-empty half-open invariant as `Interval`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IndexedInterval<T, I = u64> {
    /// Checked non-empty interval.
    pub interval: Interval<T>,
    /// External index or identifier carried alongside the interval.
    pub idx: I,
}

impl<T, I> IndexedInterval<T, I>
where
    T: Copy + PartialOrd + Display,
    I: Copy,
{
    /// Create a checked indexed interval from raw bounds and an index.
    ///
    /// This validates the interval bounds and keeps the index unchanged.
    ///
    /// Parameters
    /// ----------
    /// - `start`:
    ///   Inclusive start coordinate.
    /// - `end`:
    ///   Exclusive end coordinate. Must be greater than `start`.
    /// - `idx`:
    ///   External index or identifier to carry with the interval.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Checked indexed interval.
    pub fn new(start: T, end: T, idx: I) -> Result<Self> {
        Ok(Self {
            interval: Interval::new(start, end)?,
            idx,
        })
    }

    /// Create an indexed interval from an already checked interval.
    ///
    /// Parameters
    /// ----------
    /// - `interval`:
    ///   Checked non-empty interval.
    /// - `idx`:
    ///   External index or identifier to carry with the interval.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Indexed interval with the provided metadata.
    pub fn from_interval(interval: Interval<T>, idx: I) -> Self {
        Self { interval, idx }
    }

    /// Return the inclusive start coordinate.
    #[inline]
    pub fn start(&self) -> T {
        self.interval.start()
    }

    /// Return the exclusive end coordinate.
    #[inline]
    pub fn end(&self) -> T {
        self.interval.end()
    }

    /// Return the carried index.
    #[inline]
    pub fn idx(&self) -> I {
        self.idx
    }

    /// Return the interval and index as `(start, end, idx)`.
    #[inline]
    pub fn into_tuple(self) -> (T, T, I) {
        (self.interval.start(), self.interval.end(), self.idx)
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

impl<T, I> TryFrom<(T, T, I)> for IndexedInterval<T, I>
where
    T: Copy + PartialOrd + Display,
    I: Copy,
{
    type Error = Error;

    /// Convert a `(start, end, idx)` tuple into a checked indexed interval.
    ///
    /// Use this when interval bounds and their external identifier already
    /// exist as tuples and should be validated during conversion.
    ///
    /// Parameters
    /// ----------
    /// - `bounds`:
    ///   Interval bounds and index as `(start, end, idx)`.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Checked indexed interval.
    ///
    /// Example
    /// -------
    /// ```rust
    /// use cfdnalab::shared::interval::IndexedInterval;
    ///
    /// let windows: cfdnalab::Result<Vec<_>> = vec![(5_u64, 6_u64, 10_u64), (10, 20, 11)]
    ///     .into_iter()
    ///     .map(IndexedInterval::try_from)
    ///     .collect();
    ///
    /// assert_eq!(windows?[0].idx(), 10);
    /// # Ok::<(), cfdnalab::Error>(())
    /// ```
    fn try_from(bounds: (T, T, I)) -> Result<Self> {
        Self::new(bounds.0, bounds.1, bounds.2)
    }
}
