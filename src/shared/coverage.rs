use crate::shared::fragment::{minimal_fragment::Fragment, segment_fragment::FragmentWithSegments};
use crate::shared::interval::Interval;
use anyhow::{Result, bail};
use rayon::prelude::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Stage {
    Building, // Accepting +w/-w edits in delta
    Covered,  // coverage present, indexes may or may not be built
    Indexed,  // coverage present, indexes built
}

/// Dense per-base coverage with optional blacklist and O(1) interval queries via prefix sums.
///
/// This collects fragments into a +w/-w delta array, converts the delta to per-base coverage,
/// and builds prefix-sum indexes so you can query sums and averages over any interval.
///
/// Example
/// -------
/// ```rust
/// use cfdnalab::shared::coverage::Coverage;
/// use cfdnalab::shared::fragment::minimal_fragment::Fragment;
/// use cfdnalab::shared::gc_tag::GCTagValue;
/// use cfdnalab::shared::interval::Interval;
///
/// # use anyhow::Result;
/// # fn demo() -> Result<()> {
/// let length: u32 = 1_000_000; // e.g., chrom_len
/// let mut cp = Coverage::new(length);
///
/// // Unweighted fragment
/// cp.add_fragment(Fragment {
///     tid: 0,
///     interval: Interval::new(100, 200)?,
///     gc_tag: GCTagValue::default(),
/// })?;
///
/// // GC-weighted fragment
/// cp.add_fragment_weighted(
///     Fragment {
///         tid: 0,
///         interval: Interval::new(150, 250)?,
///         gc_tag: GCTagValue::default(),
///     },
///     0.87,
/// )?;
///
/// // Optional blacklist
/// let blacklist_intervals = Interval::from_tuples(&[(120, 140), (150, 153)])?;
/// cp.set_blacklist_mask(&blacklist_intervals)?;
///
/// // Build per-base coverage and query indexes
/// cp.finalize_coverage(true);  // free delta after building coverage
/// cp.build_indexes(true)?; // build psums and free coverage
///
/// // Query averages
/// let avg_all = cp.avg_coverage(Interval::new(100, 300)?, false)?;   // Includes blacklisted bases
/// let avg_ok  = cp.avg_coverage(Interval::new(100, 300)?, true)?;    // Excludes blacklisted bases
///
/// // Raw positional coverage if needed
/// let cov = cp.coverage().unwrap();
/// assert_eq!(cov.len() as u32, length);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct Coverage {
    length: u32,                // Total sequence length in bases (e.g., chrom_len)
    delta: Vec<f64>,            // +w at start, -w at end, length = length + 1 (last is sentinel)
    coverage: Option<Vec<f32>>, // Per-base coverage after finalize_coverage, length = length
    bl_mask: Option<Vec<u8>>, // Per-base blacklist mask after finalize_blacklist_prefix, 1 = blacklisted

    // Prefix sums for fast queries
    psum_all: Option<Vec<f64>>,            // Σ coverage
    psum_unmasked: Option<Vec<f64>>,       // Σ coverage over non-blacklisted positions
    psum_unmasked_count: Option<Vec<u32>>, // Σ 1 over non-blacklisted positions

    cov_stage: Stage, // Lifecycle for coverage
}

/// Clamp finite coverage values below `floor` to exact zero.
///
/// Coverage is semantically nonnegative, so this also clamps any negative
/// floating-point artefacts to `0.0`.
pub(crate) fn clamp_finite_coverage_below_to_zero(values: &mut [f32], floor: f32) {
    let floor = floor as f64;
    for value in values.iter_mut() {
        if value.is_finite() && (*value as f64) < floor {
            *value = 0.0;
        }
    }
}

impl Coverage {
    /// Create new `Coverage` instance.
    ///
    /// Parameters
    /// ----------
    /// - length:
    ///     Total sequence length in bases (e.g., chrom_len).
    ///
    /// Returns
    /// -------
    /// - self:
    ///     New empty prefix.
    pub fn new(length: u32) -> Self {
        Self {
            length,
            delta: vec![0.0_f64; length as usize + 1],
            coverage: None,
            bl_mask: None,
            psum_all: None,
            psum_unmasked: None,
            psum_unmasked_count: None,
            cov_stage: Stage::Building,
        }
    }

    /// add_fragment: +1 at start, -1 at end
    ///
    /// Parameters
    /// ----------
    /// - frag:
    ///     Fragment on the reference `[start, end)`, 0-based, end-exclusive.
    #[inline]
    pub fn add_fragment(&mut self, frag: Fragment) -> Result<()> {
        self.add_fragment_weighted(frag, 1.0)
    }

    /// add_fragment with floating weight w: +w at start, -w at end
    ///
    /// Parameters
    /// ----------
    /// - frag:
    ///     Fragment on the reference `[start, end)`.
    /// - weight:
    ///     Weight to add, must be finite and >= 0.
    #[inline]
    pub fn add_fragment_weighted(&mut self, frag: Fragment, weight: f64) -> Result<()> {
        if !self.prefix_available() {
            anyhow::bail!(
                "prefix was dropped; cannot add fragments. Rebuild or create a new Coverage"
            );
        }
        if !weight.is_finite() || weight < 0.0 {
            anyhow::bail!("invalid weight {}", weight);
        }

        // If we are finalized, invalidate coverage and indexes and go back to Building
        if matches!(self.cov_stage, Stage::Covered | Stage::Indexed) {
            self.coverage = None;
            self.invalidate_indexes();
            self.cov_stage = Stage::Building;
        }

        let n = self.delta.len();
        let start = frag.start() as usize;
        let end = frag.end() as usize;

        if end > self.length as usize || end >= n {
            anyhow::bail!(
                "fragment end {} out of bounds for sequence length {}",
                end,
                self.length
            );
        }

        // Apply boundary deltas
        self.delta[start] += weight;
        self.delta[end] -= weight;

        // Invalidate indexes because the underlying data changed
        self.invalidate_indexes();
        Ok(())
    }

    /// add_fragment_with_segments: add a fragment using either its full span or explicit segments
    ///
    /// Summary
    /// -------
    /// Accepts a `FragmentWithSegments` that already encodes any desired behavior:
    /// - If `segments` is `None`, we add the plain fragment span `[start, end)`
    /// - If `segments` is `Some`, we add those `[start, end)` segments (already unioned and
    ///   optionally including the inter-mate gap if requested upstream)
    ///
    /// Notes
    /// -----
    /// Inter-mate gap handling and ref-gap segmentation are decided upstream in
    /// `collect_fragment_with_segments`, so this method only applies what it is given.
    pub fn add_fragment_with_segments(
        &mut self,
        frag: FragmentWithSegments,
        weight: f64,
    ) -> anyhow::Result<()> {
        if !self.prefix_available() {
            anyhow::bail!(
                "prefix was dropped; cannot add fragments. Rebuild or create a new Coverage"
            );
        }
        if !weight.is_finite() || weight < 0.0 {
            anyhow::bail!("invalid weight {}", weight);
        }
        // If coverage/indexes exist, invalidate and go back to Building
        if matches!(self.cov_stage, Stage::Covered | Stage::Indexed) {
            self.coverage = None;
            self.invalidate_indexes();
            self.cov_stage = Stage::Building;
        }

        match frag.segments {
            None => {
                // Plain span
                let base = Fragment {
                    tid: frag.tid,
                    interval: frag.interval,
                    gc_tag: crate::shared::gc_tag::GCTagValue::default(),
                };
                self.add_fragment_weighted(base, weight)
            }
            Some(segs) => {
                // Apply +w/-w per segment
                let n = self.delta.len();
                let len = self.length as usize;
                for segment in segs {
                    let start_idx = segment.start() as usize;
                    let end_idx = segment.end() as usize;
                    if end_idx > len || start_idx >= n {
                        anyhow::bail!(
                            "segment [{}..{}) out of bounds for sequence length {}",
                            segment.start(),
                            segment.end(),
                            self.length
                        );
                    }
                    self.delta[start_idx] += weight;
                    self.delta[end_idx] -= weight;
                }
                self.invalidate_indexes();
                Ok(())
            }
        }
    }

    /// finalize_coverage: build per-base coverage from the +w/-w prefix.
    ///
    /// If `drop_delta` is true, the +w/-w `delta` is freed at the end.
    ///
    /// Returns
    /// -------
    /// - coverage:
    ///     Borrowed slice of per-base coverage with length = `length`.
    pub fn finalize_coverage(&mut self, drop_delta: bool) -> &[f32] {
        // Build coverage from the +w/-w prefix without destroying delta
        let mut cov = vec![0.0_f32; self.length as usize];

        // Cumulative sum over delta
        let mut run = 0.0_f64;
        for i in 0..=self.length as usize {
            run += self.delta[i];
            if i < self.length as usize {
                cov[i] = run as f32;
            }
        }

        self.coverage = Some(cov);
        self.invalidate_indexes();
        self.cov_stage = Stage::Covered;

        if drop_delta {
            self.drop_deltas();
        }
        self.coverage.as_ref().unwrap()
    }

    /// Set or replace the blacklist mask from half-open intervals `[start, end)`,
    /// expressed in the **same coordinate space as this `Coverage`**
    /// (i.e., prefix-local `0..self.length`).
    ///
    /// Contract
    /// - `start < end`
    /// - `0 <= start` and `end <= self.length`
    /// - Intervals may overlap; overlaps are allowed and merged.
    /// - If `intervals` is empty, the blacklist mask is removed (`None`) to avoid
    ///   allocating an all-zero vector.
    ///
    /// Errors
    /// - Returns an error if any interval violates the contract (out of bounds or empty).
    pub fn set_blacklist_mask(&mut self, intervals: &[Interval<u64>]) -> Result<()> {
        if intervals.is_empty() {
            // No blacklist -> drop mask to avoid allocating an all-zero vector.
            self.bl_mask = None;
            self.invalidate_indexes();
            return Ok(());
        }

        let n = self.length as usize;
        let mut mask = vec![0u8; n];

        for interval in intervals {
            let s64 = interval.start();
            let e64 = interval.end();
            if s64 >= e64 {
                bail!("blacklist interval start {} >= end {}", s64, e64);
            }
            if e64 > self.length as u64 {
                bail!(
                    "out of bounds: blacklist interval end {} exceeds sequence length {}",
                    e64,
                    self.length
                );
            }
            // Safe after the checks above
            let a = s64 as usize;
            let b = e64 as usize;
            debug_assert!(b <= n);
            mask[a..b].fill(1);
        }

        self.bl_mask = Some(mask);
        self.invalidate_indexes(); // sums depend on mask
        Ok(())
    }

    /// build_indexes: prepare prefix sums for fast interval queries
    ///
    /// What gets built
    /// ---------------
    /// - Always builds `psum_all`  (Σ coverage over all bases) with length = n+1
    /// - Only builds `psum_unmasked` and `psum_unmasked_count` when a blacklist mask exists
    ///
    /// Safety
    /// ------
    /// - Requires `finalize_coverage()` to have been called
    pub fn build_indexes(&mut self, drop_coverage: bool) -> anyhow::Result<()> {
        let cov = self.coverage.as_ref().ok_or_else(|| {
            anyhow::anyhow!("coverage not finalized, call finalize_coverage() first")
        })?;
        let n = cov.len();

        // Always build psum_all (n+1 with empty prefix at index 0)
        let mut psum_all = Vec::with_capacity(n + 1);
        psum_all.push(0.0_f64);

        if let Some(mask) = self.bl_mask.as_ref() {
            anyhow::ensure!(
                mask.len() == n,
                "mask length {} != coverage length {}",
                mask.len(),
                n
            );
            // Mask present -> also build unmasked sums & counts
            let mut psum_unmasked = Vec::with_capacity(n + 1);
            let mut psum_unmasked_count = Vec::with_capacity(n + 1);
            psum_unmasked.push(0.0_f64);
            psum_unmasked_count.push(0u32);

            let (mut all, mut unmasked_sum, mut unmasked_count) = (0.0_f64, 0.0_f64, 0u32);
            for i in 0..n {
                let c = cov[i] as f64;
                all += c;
                psum_all.push(all);

                if mask[i] == 0 {
                    unmasked_sum += c;
                    unmasked_count = unmasked_count.saturating_add(1);
                }
                psum_unmasked.push(unmasked_sum);
                psum_unmasked_count.push(unmasked_count);
            }

            self.psum_unmasked = Some(psum_unmasked);
            self.psum_unmasked_count = Some(psum_unmasked_count);
        } else {
            // No mask -> only psum_all; keep unmasked structures None to save RAM
            let mut all = 0.0_f64;
            for i in 0..n {
                all += cov[i] as f64;
                psum_all.push(all);
            }
            self.psum_unmasked = None;
            self.psum_unmasked_count = None;
        }

        self.psum_all = Some(psum_all);
        self.cov_stage = Stage::Indexed;
        if drop_coverage {
            self.drop_coverage();
        }
        Ok(())
    }

    /// Return the coverage sum over an interval.
    ///
    /// Parameters
    /// ----------
    /// - `interval`:
    ///     Checked half-open interval `[start, end)` to query.
    /// - `exclude_blacklisted`:
    ///     Exclude blacklisted positions from the sum when a blacklist mask is available.
    ///
    /// Returns
    /// -------
    /// - `sum`:
    ///     Coverage sum over the interval, masked if requested.
    ///
    /// Notes
    /// -----
    /// - `interval` must be non-empty because `Interval` enforces `end > start`.
    /// - The interval must lie within `0..self.length`.
    #[inline]
    pub fn sum_coverage(
        &mut self,
        interval: Interval<u32>,
        exclude_blacklisted: bool,
    ) -> Result<f64> {
        self.ensure_ready_for_queries()?;
        self.check_interval(interval)?;

        let start_idx = interval.start() as usize;
        let end_idx = interval.end() as usize;

        let sum = if exclude_blacklisted && self.has_blacklist() {
            let prefix_sums_unmasked = self.psum_unmasked.as_ref().unwrap();
            prefix_sums_unmasked[end_idx] - prefix_sums_unmasked[start_idx]
        } else {
            let prefix_sums_all = self.psum_all.as_ref().unwrap();
            prefix_sums_all[end_idx] - prefix_sums_all[start_idx]
        };
        Ok(sum)
    }

    /// Return the average coverage over an interval.
    ///
    /// Parameters
    /// ----------
    /// - `interval`:
    ///     Checked half-open interval `[start, end)` to query.
    /// - `exclude_blacklisted`:
    ///     Exclude blacklisted positions from the average when a blacklist mask is available.
    ///
    /// Returns
    /// -------
    /// - `avg`:
    ///     Coverage average over the interval, masked if requested.
    ///
    /// Notes
    /// -----
    /// - `interval` must be non-empty because `Interval` enforces `end > start`.
    /// - When `exclude_blacklisted` is true and every position in the interval is blacklisted,
    ///   this returns `0.0`.
    #[inline]
    pub fn avg_coverage(
        &mut self,
        interval: Interval<u32>,
        exclude_blacklisted: bool,
    ) -> Result<f32> {
        self.ensure_ready_for_queries()?;
        self.check_interval(interval)?;

        let start_idx = interval.start() as usize;
        let end_idx = interval.end() as usize;

        if exclude_blacklisted && self.has_blacklist() {
            let prefix_sums_unmasked = self.psum_unmasked.as_ref().unwrap();
            let prefix_unmasked_counts = self.psum_unmasked_count.as_ref().unwrap();
            let sum = prefix_sums_unmasked[end_idx] - prefix_sums_unmasked[start_idx];
            let unmasked_position_count =
                prefix_unmasked_counts[end_idx] - prefix_unmasked_counts[start_idx];
            if unmasked_position_count == 0 {
                return Ok(0.0);
            }
            Ok((sum / unmasked_position_count as f64) as f32)
        } else {
            let prefix_sums_all = self.psum_all.as_ref().unwrap();
            let sum = prefix_sums_all[end_idx] - prefix_sums_all[start_idx];
            Ok((sum / interval.len() as f64) as f32)
        }
    }

    /// Return coverage sums for many intervals using prefix sums.
    ///
    /// Parameters
    /// ----------
    /// - `intervals`:
    ///     Checked half-open intervals `[start, end)` to query.
    /// - `exclude_blacklisted`:
    ///     Exclude blacklisted positions from the sums when a blacklist mask is available.
    /// - `parallelize`:
    ///     Process intervals with rayon parallel iterators.
    ///
    /// Returns
    /// -------
    /// - `sums`:
    ///     One coverage sum per interval, in the same order as the input slice.
    #[inline]
    pub fn bulk_sum_coverage(
        &mut self,
        intervals: &[Interval<u32>],
        exclude_blacklisted: bool,
        parallelize: bool,
    ) -> Result<Vec<f64>> {
        self.ensure_ready_for_queries()?;
        for &interval in intervals {
            self.check_interval(interval)?;
        }

        Ok(if exclude_blacklisted && self.has_blacklist() {
            let prefix_sums_unmasked = self.psum_unmasked.as_ref().unwrap();
            self.bulk_sums_with_prefix_intervals(intervals, prefix_sums_unmasked, parallelize)
        } else {
            let prefix_sums_all = self.psum_all.as_ref().unwrap();
            self.bulk_sums_with_prefix_intervals(intervals, prefix_sums_all, parallelize)
        })
    }

    /// Return coverage averages for many intervals using prefix sums.
    ///
    /// Parameters
    /// ----------
    /// - `intervals`:
    ///     Checked half-open intervals `[start, end)` to query.
    /// - `exclude_blacklisted`:
    ///     Exclude blacklisted positions from the averages when a blacklist mask is available.
    /// - `parallelize`:
    ///     Process intervals with rayon parallel iterators.
    ///
    /// Returns
    /// -------
    /// - `avgs`:
    ///     One average coverage value per interval, in the same order as the input slice.
    ///
    /// Notes
    /// -----
    /// - When `exclude_blacklisted` is true and an interval is fully blacklisted, the returned
    ///   average for that interval is `0.0`.
    #[inline]
    pub fn bulk_avg_coverage(
        &mut self,
        intervals: &[Interval<u32>],
        exclude_blacklisted: bool,
        parallelize: bool,
    ) -> Result<Vec<f32>> {
        self.ensure_ready_for_queries()?;
        for &interval in intervals {
            self.check_interval(interval)?;
        }

        Ok(if exclude_blacklisted && self.has_blacklist() {
            let prefix_sums_unmasked = self.psum_unmasked.as_ref().unwrap();
            let prefix_unmasked_counts = self.psum_unmasked_count.as_ref().unwrap();
            self.bulk_avgs_with_prefix_intervals(
                intervals,
                prefix_sums_unmasked,
                Some(prefix_unmasked_counts),
                parallelize,
            )
        } else {
            let prefix_sums_all = self.psum_all.as_ref().unwrap();
            self.bulk_avgs_with_prefix_intervals(intervals, prefix_sums_all, None, parallelize)
        })
    }

    /// Return raw coverage at the requested positions.
    ///
    /// Parameters
    /// ----------
    /// - positions:
    ///     Positions to fetch, 0-based.
    ///
    /// Returns
    /// -------
    /// - values:
    ///     Coverage values at each requested position.
    pub fn coverage_at_positions(&self, positions: &[u32]) -> Result<Vec<f32>> {
        let cov = match self.coverage.as_ref() {
            Some(c) => c,
            None => {
                anyhow::bail!("coverage not finalized, call finalize_coverage() first")
            }
        };

        let len = self.length as usize;
        let mut out = Vec::with_capacity(positions.len());

        for &p in positions {
            let i = p as usize;
            if i >= len {
                anyhow::bail!("position {} out of bounds for length {}", p, self.length);
            }
            out.push(cov[i]);
        }
        Ok(out)
    }

    /// Return coverage at positions, with blacklisted sites as `NaN`.
    ///
    /// When no blacklist mask is present, no elements will be `NaN`.
    ///
    /// Parameters
    /// ----------
    /// - positions:
    ///     Positions to fetch, 0-based.
    ///
    /// Returns
    /// -------
    /// - values:
    ///     Coverage values; blacklisted sites are `f32::NAN`.
    pub fn coverage_at_positions_nan(&self, positions: &[u32]) -> Result<Vec<f32>> {
        let cov = match self.coverage.as_ref() {
            Some(c) => c,
            None => {
                anyhow::bail!("coverage not finalized, call finalize_coverage() first")
            }
        };

        let len = self.length as usize;
        let mask_opt = self.bl_mask.as_ref();
        let mut out = Vec::with_capacity(positions.len());

        match mask_opt {
            Some(mask) => {
                for &p in positions {
                    let i = p as usize;
                    if i >= len {
                        anyhow::bail!("position {} out of bounds for length {}", p, self.length);
                    }
                    out.push(if mask[i] == 0 { cov[i] } else { f32::NAN });
                }
            }
            None => {
                for &p in positions {
                    let i = p as usize;
                    if i >= len {
                        anyhow::bail!("position {} out of bounds for length {}", p, self.length);
                    }
                    out.push(cov[i]);
                }
            }
        }
        Ok(out)
    }

    /// Return blacklist mask values at the requested positions.
    ///
    /// Value is 1 for blacklisted, 0 for unmasked.
    ///
    /// Parameters
    /// ----------
    /// - positions:
    ///     Positions to fetch, 0-based.
    ///
    /// Returns
    /// -------
    /// - mask:
    ///     Mask values at each requested position; if no blacklist, all zeros.
    pub fn mask_at_positions(&self, positions: &[u32]) -> Result<Vec<u8>> {
        let len = self.length as usize;
        let mut out = Vec::with_capacity(positions.len());

        match self.bl_mask.as_ref() {
            Some(mask) => {
                for &p in positions {
                    let i = p as usize;
                    if i >= len {
                        anyhow::bail!("position {} out of bounds for length {}", p, self.length);
                    }
                    out.push(mask[i]);
                }
            }
            None => {
                for &p in positions {
                    let i = p as usize;
                    if i >= len {
                        anyhow::bail!("position {} out of bounds for length {}", p, self.length);
                    }
                    out.push(0u8);
                }
            }
        }
        Ok(out)
    }

    /// Return positional coverage for `[start, end)`
    ///
    /// Parameters
    /// ----------
    /// - start:
    ///     Start of interval, inclusive
    /// - end:
    ///     End of interval, exclusive
    /// - nan_blacklisted:
    ///     Set positions overlapped by the blacklist to `f32::NAN`
    ///
    /// Returns
    /// -------
    /// - values:
    ///     Coverage values for each position in `[start, end)`, with optional `f32::NAN`s for blacklisted
    pub fn coverage_in_window(
        &self,
        start: u32,
        end: u32,
        nan_blacklisted: bool,
    ) -> anyhow::Result<Vec<f32>> {
        self.ensure_coverage()?;

        // Bounds check
        self.check_bounds(start, end)?;

        let cov = self.coverage.as_ref().unwrap();
        let a = start as usize;
        let b = end as usize;

        if !nan_blacklisted {
            return Ok(cov[a..b].to_vec());
        }

        // Apply NaN to blacklisted positions if mask exists
        match self.bl_mask.as_ref() {
            Some(mask) => {
                let mut out = Vec::with_capacity(b - a);
                for i in a..b {
                    out.push(if mask[i] == 1 { f32::NAN } else { cov[i] });
                }
                Ok(out)
            }
            None => Ok(cov[a..b].to_vec()),
        }
    }

    /// coverage: borrowed per-base coverage slice if finalized
    ///
    /// Returns
    /// -------
    /// - coverage:
    ///     Per-base coverage of length `length`, if available.
    pub fn coverage(&self) -> Option<&[f32]> {
        self.coverage.as_deref()
    }

    // coverage: mutable borrowed per-base coverage slice if finalized
    ///
    /// Returns
    /// -------
    /// - coverage:
    ///     Mutable per-base coverage of length `length`, if available.
    pub fn coverage_mut(&mut self) -> Option<&mut [f32]> {
        self.coverage.as_deref_mut()
    }

    /// blacklist_mask: borrowed per-base mask slice if finalized
    ///
    /// Returns
    /// -------
    /// - mask:
    ///     Per-base mask where 1 = blacklisted.
    pub fn blacklist_mask(&self) -> Option<&[u8]> {
        self.bl_mask.as_deref()
    }

    /// Length accessor
    pub fn length(&self) -> u32 {
        self.length
    }

    /// Length acessor (alias for length)
    pub fn len(&self) -> u32 {
        self.length()
    }

    /// Drop the +w/-w delta to free memory. Further add_* calls will error.
    pub fn drop_deltas(&mut self) {
        self.delta.clear();
        self.delta.shrink_to_fit();
    }

    /// Drop the coverage vector to free memory. The various coverage getters will error.
    pub fn drop_coverage(&mut self) {
        self.coverage = None;
    }

    // Remove the current blacklist if it exists
    pub fn clear_blacklist(&mut self) {
        self.bl_mask = None;
        self.invalidate_indexes();
    }

    #[inline]
    pub fn psum_all_ref(&self) -> Option<&[f64]> {
        self.psum_all.as_deref()
    }
    #[inline]
    pub fn psum_unmasked_ref(&self) -> Option<&[f64]> {
        self.psum_unmasked.as_deref()
    }
    #[inline]
    pub fn psum_unmasked_count_ref(&self) -> Option<&[u32]> {
        self.psum_unmasked_count.as_deref()
    }

    // Private helpers

    #[inline]
    fn ensure_ready_for_queries(&mut self) -> anyhow::Result<()> {
        match self.cov_stage {
            Stage::Indexed => Ok(()),                    // good to go
            Stage::Covered => self.build_indexes(false), // lazily build, then OK
            Stage::Building => {
                anyhow::bail!("coverage not finalized; call finalize_coverage() first")
            }
        }
    }

    fn invalidate_indexes(&mut self) {
        self.psum_all = None;
        self.psum_unmasked = None;
        self.psum_unmasked_count = None;
    }

    #[inline]
    fn ensure_coverage(&self) -> Result<()> {
        if self.coverage.is_none() {
            anyhow::bail!("coverage not finalized; call finalize_coverage() first");
        }
        Ok(())
    }

    fn check_bounds(&self, start: u32, end: u32) -> Result<()> {
        if start > end {
            anyhow::bail!("start {} > end {}", start, end);
        }
        if end > self.length {
            anyhow::bail!("end {} exceeds sequence length {}", end, self.length);
        }
        Ok(())
    }

    fn check_interval(&self, interval: Interval<u32>) -> Result<()> {
        if interval.end() > self.length {
            anyhow::bail!(
                "end {} exceeds sequence length {}",
                interval.end(),
                self.length
            );
        }
        Ok(())
    }

    fn bulk_sums_with_prefix_intervals(
        &self,
        intervals: &[Interval<u32>],
        prefix_sums: &[f64],
        parallelize: bool,
    ) -> Vec<f64> {
        let compute = |interval: &Interval<u32>| -> f64 {
            let start_idx = interval.start() as usize;
            let end_idx = interval.end() as usize;
            prefix_sums[end_idx] - prefix_sums[start_idx]
        };

        if parallelize {
            intervals.par_iter().map(compute).collect()
        } else {
            intervals.iter().map(compute).collect()
        }
    }

    fn bulk_avgs_with_prefix_intervals(
        &self,
        intervals: &[Interval<u32>],
        prefix_sums: &[f64],
        unmasked_counts: Option<&[u32]>,
        parallelize: bool,
    ) -> Vec<f32> {
        let compute = |interval: &Interval<u32>| -> f32 {
            let start_idx = interval.start() as usize;
            let end_idx = interval.end() as usize;
            let sum = prefix_sums[end_idx] - prefix_sums[start_idx];

            if let Some(counts) = unmasked_counts {
                let unmasked_position_count = counts[end_idx] - counts[start_idx];
                if unmasked_position_count == 0 {
                    0.0
                } else {
                    (sum / unmasked_position_count as f64) as f32
                }
            } else {
                (sum / interval.len() as f64) as f32
            }
        };

        if parallelize {
            intervals.par_iter().map(compute).collect()
        } else {
            intervals.iter().map(compute).collect()
        }
    }

    #[inline]
    fn prefix_available(&self) -> bool {
        !self.delta.is_empty()
    }

    // True if a blacklist is present
    pub fn has_blacklist(&self) -> bool {
        self.bl_mask.is_some()
    }
}

#[cfg(test)]
mod tests {
    include!("coverage_tests.rs");
    include!("test_coverage_correlation.rs");
}
