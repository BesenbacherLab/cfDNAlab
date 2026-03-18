use crate::{Error, Result};
use num_traits::{CheckedAdd, CheckedSub, Signed, Unsigned};
use std::fmt::Display;
use std::ops::Sub;

/// A checked half-open interval `[start, end)`.
///
/// Use this for the geometric part of domain structs that carry genomic spans.
/// This type only represents non-empty intervals, so construction requires
/// `end > start`. For ordered bounds that may be empty, use `Span<T>`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Interval<T> {
    start: T,
    end: T,
}

/// Controls whether interval-list merging collapses touching intervals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TouchingMergePolicy {
    /// Merge only intervals with positive overlap.
    KeepTouchingSeparate,
    /// Merge intervals that overlap or meet at one boundary.
    MergeTouching,
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
    pub fn as_tuple(&self) -> (T, T) {
        (self.start, self.end)
    }

    /// Return the interval bounds as `(start, end)`.
    #[inline]
    pub fn into_inner(self) -> (T, T) {
        (self.start, self.end)
    }

    /// Convert a slice of `(start, end)` tuples into checked intervals.
    ///
    /// Use this when helpers or fixtures already store genomic spans as tuples and
    /// you want one checked conversion step before passing them into APIs that use
    /// `Interval`.
    ///
    /// Parameters
    /// ----------
    /// - `entries`:
    ///   Slice of `(start, end)` tuples to validate and convert.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Vector of checked intervals in the same order as the input slice.
    ///
    /// Example
    /// -------
    /// ```rust
    /// use cfdnalab::shared::interval::Interval;
    ///
    /// let intervals = Interval::from_tuples(&[(5_u64, 10_u64), (10, 20)])?;
    ///
    /// assert_eq!(intervals.len(), 2);
    /// assert_eq!(intervals[1].start(), 10);
    /// # Ok::<(), cfdnalab::Error>(())
    /// ```
    pub fn from_tuples(entries: &[(T, T)]) -> Result<Vec<Self>> {
        entries
            .iter()
            .map(|&(start, end)| Self::new(start, end))
            .collect()
    }
}

impl<T> Interval<T>
where
    T: Copy + Sub<Output = T>,
{
    /// Return the interval length as `end - start`.
    ///
    /// Because this type only allows non-empty half-open intervals, the result
    /// is always greater than zero for numeric coordinate types.
    #[inline]
    pub fn len(&self) -> T {
        self.end - self.start
    }
}

impl<T> Interval<T>
where
    T: Copy + Ord,
{
    /// Return whether `point` lies inside this half-open interval.
    #[inline]
    pub fn contains_point(&self, point: T) -> bool {
        point >= self.start && point < self.end
    }

    /// Return whether `other` lies fully inside this interval.
    #[inline]
    pub fn contains_interval(&self, other: Self) -> bool {
        other.start >= self.start && other.end <= self.end
    }

    /// Return whether two half-open intervals intersect.
    #[inline]
    pub fn intersects(&self, other: Self) -> bool {
        other.end > self.start && other.start < self.end
    }

    /// Return a new interval holding the shared part of two half-open intervals, if any.
    #[inline]
    pub fn intersection(&self, other: Self) -> Option<Self> {
        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        (end > start).then_some(Self { start, end })
    }

    /// Return a new interval clipped to `bounds`, if any span remains.
    ///
    /// Use this when the receiver is the primary interval and `bounds` provides
    /// the allowed coordinate range. This is equivalent to the interval
    /// intersection, but reads more naturally at call sites that are clipping a
    /// value to enclosing bounds. This does not mutate the receiver.
    #[inline]
    pub fn clip_to(&self, bounds: Self) -> Option<Self> {
        self.intersection(bounds)
    }

    /// Return a new interval clipped so it starts no earlier than `lower_bound`.
    ///
    /// This does not mutate the receiver. Returns `None` when the clipped
    /// interval would be empty.
    #[inline]
    pub fn clip_lower(&self, lower_bound: T) -> Option<Self> {
        let start = self.start.max(lower_bound);
        (self.end > start).then_some(Self {
            start,
            end: self.end,
        })
    }

    /// Return a new interval clipped so it ends no later than `upper_bound`.
    ///
    /// This does not mutate the receiver. Returns `None` when the clipped
    /// interval would be empty.
    #[inline]
    pub fn clip_upper(&self, upper_bound: T) -> Option<Self> {
        let end = self.end.min(upper_bound);
        (end > self.start).then_some(Self {
            start: self.start,
            end,
        })
    }

    /// Return a new interval spanning both inputs.
    ///
    /// This does not mutate either operand. Because both operands are already
    /// checked non-empty intervals, the spanning interval is also guaranteed to
    /// be a valid non-empty half-open interval.
    #[inline]
    pub fn expand_to_include(&self, other: Self) -> Self {
        Self {
            start: self.start.min(other.start),
            end: self.end.max(other.end),
        }
    }
}

/// Append `next_interval` to a start-sorted interval list, merging with the previous
/// entry when the chosen touching policy allows it.
///
/// Use this when you are building a merged interval list incrementally and the
/// input order is already sorted by interval start. The merge policy controls
/// whether touching intervals like `[10,20)` and `[20,30)` collapse into one
/// interval or remain separate.
pub fn push_merged_interval<T>(
    merged_intervals: &mut Vec<Interval<T>>,
    next_interval: Interval<T>,
    touching_merge_policy: TouchingMergePolicy,
) where
    T: Copy + Ord,
{
    if let Some(last_interval) = merged_intervals.last_mut() {
        let should_merge = match touching_merge_policy {
            TouchingMergePolicy::KeepTouchingSeparate => last_interval.end > next_interval.start,
            TouchingMergePolicy::MergeTouching => last_interval.end >= next_interval.start,
        };
        if should_merge {
            if next_interval.end > last_interval.end {
                *last_interval = last_interval.expand_to_include(next_interval);
            }
            return;
        }
    }

    merged_intervals.push(next_interval);
}

/// Merge a start-sorted interval list with the chosen touching policy.
///
/// This assumes `intervals` are already sorted by `(start, end)`. Use this when
/// you already have the full interval vector and want the merged result without
/// re-sorting it.
pub fn merge_sorted_intervals<T>(
    intervals: Vec<Interval<T>>,
    touching_merge_policy: TouchingMergePolicy,
) -> Vec<Interval<T>>
where
    T: Copy + Ord,
{
    let mut merged_intervals: Vec<Interval<T>> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        push_merged_interval(&mut merged_intervals, interval, touching_merge_policy);
    }
    merged_intervals
}

/// Sort and merge an interval list with the chosen touching policy.
///
/// Use this when you already have the full interval vector but do not know
/// whether it is sorted by `(start, end)`.
pub fn merge_intervals<T>(
    mut intervals: Vec<Interval<T>>,
    touching_merge_policy: TouchingMergePolicy,
) -> Vec<Interval<T>>
where
    T: Copy + Ord,
{
    intervals.sort_unstable_by_key(|interval| (interval.start, interval.end));
    merge_sorted_intervals(intervals, touching_merge_policy)
}

impl<T> Interval<T>
where
    T: Copy + Display + CheckedAdd + Signed,
{
    /// Return a new interval shifted by `delta`.
    ///
    /// The offset uses the same signed numeric type as the interval
    /// coordinates. A negative `delta` shifts the interval left. This does not
    /// mutate the receiver and returns an error when the shifted bounds would
    /// overflow or underflow the coordinate type.
    pub fn offset(self, delta: T) -> Result<Self> {
        let start = match self.start.checked_add(&delta) {
            Some(value) => value,
            None => {
                return Err(Error::InvalidIntervalOffset {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    offset: delta.to_string(),
                });
            }
        };
        let end = match self.end.checked_add(&delta) {
            Some(value) => value,
            None => {
                return Err(Error::InvalidIntervalOffset {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    offset: delta.to_string(),
                });
            }
        };
        Ok(Self { start, end })
    }
}

impl<T> Interval<T>
where
    T: Copy + Display + CheckedAdd + Unsigned,
{
    /// Return a new interval shifted right by `amount`.
    ///
    /// This does not mutate the receiver and returns an error when the shifted
    /// bounds would overflow the coordinate type.
    pub fn shift_right(self, amount: T) -> Result<Self> {
        let start = match self.start.checked_add(&amount) {
            Some(value) => value,
            None => {
                return Err(Error::InvalidIntervalOffset {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    offset: amount.to_string(),
                });
            }
        };
        let end = match self.end.checked_add(&amount) {
            Some(value) => value,
            None => {
                return Err(Error::InvalidIntervalOffset {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    offset: amount.to_string(),
                });
            }
        };
        Ok(Self { start, end })
    }
}

impl<T> Interval<T>
where
    T: Copy + Display + CheckedSub + Unsigned,
{
    /// Return a new interval shifted left by `amount`.
    ///
    /// This does not mutate the receiver and returns an error when the shifted
    /// bounds would underflow the coordinate type.
    pub fn shift_left(self, amount: T) -> Result<Self> {
        let start = match self.start.checked_sub(&amount) {
            Some(value) => value,
            None => {
                return Err(Error::InvalidIntervalOffset {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    offset: amount.to_string(),
                });
            }
        };
        let end = match self.end.checked_sub(&amount) {
            Some(value) => value,
            None => {
                return Err(Error::InvalidIntervalOffset {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    offset: amount.to_string(),
                });
            }
        };
        Ok(Self { start, end })
    }
}

impl Interval<u64> {
    /// Convert a checked unsigned interval into a checked signed interval.
    ///
    /// Use this when an external API requires `i64` coordinates even though the interval is
    /// represented internally with non-negative genomic coordinates.
    pub fn try_to_i64(self) -> Result<Interval<i64>> {
        let start = match i64::try_from(self.start) {
            Ok(value) => value,
            Err(_) => {
                return Err(Error::InvalidIntervalConversion {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    target_type: "i64",
                });
            }
        };
        let end = match i64::try_from(self.end) {
            Ok(value) => value,
            Err(_) => {
                return Err(Error::InvalidIntervalConversion {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    target_type: "i64",
                });
            }
        };
        Interval::new(start, end)
    }
}

impl Interval<u32> {
    /// Convert a checked unsigned interval into a checked signed interval.
    pub fn try_to_i64(self) -> Result<Interval<i64>> {
        Interval::new(self.start as i64, self.end as i64)
    }

    /// Convert a checked `u32` interval into a checked `u64` interval.
    pub fn try_to_u64(self) -> Result<Interval<u64>> {
        Interval::new(self.start as u64, self.end as u64)
    }
}

impl Interval<i64> {
    /// Convert a checked signed interval into a checked unsigned interval.
    ///
    /// Use this after signed interval arithmetic when the result must be passed
    /// back into APIs that use non-negative genomic coordinates.
    pub fn try_to_u64(self) -> Result<Interval<u64>> {
        let start = match u64::try_from(self.start) {
            Ok(value) => value,
            Err(_) => {
                return Err(Error::InvalidIntervalConversion {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    target_type: "u64",
                });
            }
        };
        let end = match u64::try_from(self.end) {
            Ok(value) => value,
            Err(_) => {
                return Err(Error::InvalidIntervalConversion {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    target_type: "u64",
                });
            }
        };
        Interval::new(start, end)
    }

    /// Convert a checked signed interval into a checked `u32` interval.
    pub fn try_to_u32(self) -> Result<Interval<u32>> {
        let start = match u32::try_from(self.start) {
            Ok(value) => value,
            Err(_) => {
                return Err(Error::InvalidIntervalConversion {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    target_type: "u32",
                });
            }
        };
        let end = match u32::try_from(self.end) {
            Ok(value) => value,
            Err(_) => {
                return Err(Error::InvalidIntervalConversion {
                    start: self.start.to_string(),
                    end: self.end.to_string(),
                    target_type: "u32",
                });
            }
        };
        Interval::new(start, end)
    }
}

/// A checked half-open span `[start, end)` that may be empty.
///
/// Use this when code needs ordered genomic bounds but the empty case is valid,
/// for example collection envelopes or optional extents. For guaranteed
/// non-empty geometry, use `Interval<T>` instead.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Span<T> {
    start: T,
    end: T,
}

impl<T> Span<T>
where
    T: Copy + PartialOrd + Display,
{
    /// Create a checked half-open span `[start, end)`.
    ///
    /// Unlike `Interval`, this type allows `end == start` for empty spans, but
    /// still rejects inverted bounds where `end < start`.
    pub fn new(start: T, end: T) -> Result<Self> {
        if end < start {
            return Err(Error::InvalidSpanBounds {
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

    /// Return the span bounds as `(start, end)`.
    #[inline]
    pub fn as_tuple(&self) -> (T, T) {
        (self.start, self.end)
    }

    /// Return the span bounds as `(start, end)`.
    #[inline]
    pub fn into_inner(self) -> (T, T) {
        (self.start, self.end)
    }

    /// Create a span from bounds that are already known to be ordered.
    #[inline]
    pub(crate) fn from_ordered(start: T, end: T) -> Self {
        debug_assert!(end >= start, "span bounds must be ordered");
        Self { start, end }
    }
}

impl<T> Span<T>
where
    T: Copy + PartialEq,
{
    /// Return whether this span is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}

impl<T> TryFrom<(T, T)> for Span<T>
where
    T: Copy + PartialOrd + Display,
{
    type Error = Error;

    /// Convert a `(start, end)` tuple into a checked half-open span.
    fn try_from(bounds: (T, T)) -> Result<Self> {
        Self::new(bounds.0, bounds.1)
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

    /// Convert a slice of `(start, end, idx)` tuples into checked indexed intervals.
    ///
    /// Use this when helpers or fixtures already store genomic windows as tuples and
    /// you want one checked conversion step before passing them into APIs that use
    /// `IndexedInterval`.
    ///
    /// Parameters
    /// ----------
    /// - `entries`:
    ///   Slice of `(start, end, idx)` tuples to validate and convert.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Vector of checked indexed intervals in the same order as the input slice.
    ///
    /// Example
    /// -------
    /// ```rust
    /// use cfdnalab::shared::interval::IndexedInterval;
    ///
    /// let windows = IndexedInterval::from_tuples(&[(5_u64, 10_u64, 0_u64), (10, 20, 1)])?;
    ///
    /// assert_eq!(windows.len(), 2);
    /// assert_eq!(windows[1].idx(), 1);
    /// # Ok::<(), cfdnalab::Error>(())
    /// ```
    pub fn from_tuples(entries: &[(T, T, I)]) -> Result<Vec<Self>> {
        entries
            .iter()
            .map(|&(start, end, idx)| Self::new(start, end, idx))
            .collect()
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
    pub fn as_tuple(&self) -> (T, T, I) {
        (self.interval.start(), self.interval.end(), self.idx)
    }

    /// Return the interval and index as `(start, end, idx)`.
    #[inline]
    pub fn into_tuple(self) -> (T, T, I) {
        (self.interval.start(), self.interval.end(), self.idx)
    }
}

impl<T, I> IndexedInterval<T, I>
where
    T: Copy + Sub<Output = T>,
    I: Copy,
{
    /// Return the interval length as `end - start`.
    ///
    /// This forwards to the checked inner interval, so the same non-empty
    /// half-open invariant applies here.
    #[inline]
    pub fn len(&self) -> T {
        self.interval.len()
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

/// A checked indexed interval together with a score or weight.
///
/// Use this when a genomic interval needs both stable caller metadata and one
/// extra numeric value, such as a score parsed from a BED-like file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScoredInterval<T, I = u64, S = f64> {
    /// Checked interval with an external index.
    pub window: IndexedInterval<T, I>,
    /// Caller-provided score associated with the interval.
    pub score: S,
}

impl<T, I, S> ScoredInterval<T, I, S>
where
    T: Copy + PartialOrd + Display,
    I: Copy,
    S: Copy,
{
    /// Create a checked scored interval from raw bounds, an index, and a score.
    pub fn new(start: T, end: T, idx: I, score: S) -> Result<Self> {
        Ok(Self {
            window: IndexedInterval::new(start, end, idx)?,
            score,
        })
    }

    /// Create a scored interval from an already checked indexed interval.
    pub fn from_indexed_interval(window: IndexedInterval<T, I>, score: S) -> Self {
        Self { window, score }
    }

    /// Convert a slice of `(start, end, idx, score)` tuples into checked scored intervals.
    pub fn from_tuples(entries: &[(T, T, I, S)]) -> Result<Vec<Self>> {
        entries
            .iter()
            .map(|&(start, end, idx, score)| Self::new(start, end, idx, score))
            .collect()
    }

    /// Return the inclusive start coordinate.
    #[inline]
    pub fn start(&self) -> T {
        self.window.start()
    }

    /// Return the exclusive end coordinate.
    #[inline]
    pub fn end(&self) -> T {
        self.window.end()
    }

    /// Return the carried index.
    #[inline]
    pub fn idx(&self) -> I {
        self.window.idx()
    }

    /// Return the carried score.
    #[inline]
    pub fn score(&self) -> S {
        self.score
    }

    /// Return the interval, index, and score as `(start, end, idx, score)`.
    #[inline]
    pub fn into_tuple(self) -> (T, T, I, S) {
        (self.start(), self.end(), self.idx(), self.score)
    }
}
