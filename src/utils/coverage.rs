use crate::utils::fragment::Fragment;
use anyhow::Result;
use rayon::prelude::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Stage {
    Building, // Accepting +w/-w edits in delta
    Covered,  // coverage present, indexes may or may not be built
    Indexed,  // coverage present, indexes built
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum BlStage {
    Absent,    // No blacklist delta present
    Building,  // Accepting +1/-1 edits in bl_delta
    Finalized, // bl_mask present (derived from bl_delta)
}

/// Prefix-based fragment coverage for a single linear sequence with optional blacklist.
///
/// This collects fragments into a +w/-w delta array, converts the delta to per-base coverage,
/// and builds prefix-sum indexes so you can query sums and averages over any interval.
///
/// Example
/// -------
/// ```rust
/// let length: u32 = 1_000_000; // e.g., chrom_len
/// let mut cp = CoveragePrefix::initialize_coverage_prefix(length);
///
/// // Unweighted fragment
/// cp.add_fragment_to_prefix(Fragment { tid: 0, start: 100, end: 200 })?;
///
/// // GC-weighted fragment
/// cp.add_fragment_to_prefix_weighted(Fragment { tid: 0, start: 150, end: 250 }, 0.87)?;
///
/// // Optional blacklist
/// cp.initialize_blacklist_prefix();
/// cp.add_blacklist_to_prefix(120, 140)?;
/// cp.finalize_blacklist_prefix();
///
/// // Build per-base coverage and query indexes
/// cp.finalize_coverage();
/// cp.build_query_index()?;
///
/// // Query averages
/// let avg_all = cp.avg_coverage(100, 300, false)?;   // Includes blacklisted bases
/// let avg_ok  = cp.avg_coverage(100, 300, true)?;    // Excludes blacklisted bases
///
/// // Raw positional coverage if needed
/// let cov = cp.coverage().unwrap();
/// assert_eq!(cov.len() as u32, length);
/// ```
#[derive(Debug, Clone)]
pub struct CoveragePrefix {
    length: u32,                // Total sequence length in bases (e.g., chrom_len)
    delta: Vec<f32>,            // +w at start, -w at end, length = length + 1 (last is sentinel)
    bl_delta: Option<Vec<i32>>, // Optional +1/-1 delta for blacklist intervals
    coverage: Option<Vec<f32>>, // Per-base coverage after finalize_coverage, length = length
    bl_mask: Option<Vec<u8>>, // Per-base blacklist mask after finalize_blacklist_prefix, 1 = blacklisted

    // Prefix sums for fast queries
    psum_all: Option<Vec<f64>>,           // Σ coverage
    psum_allowed: Option<Vec<f64>>,       // Σ coverage over non-blacklisted positions
    psum_allowed_count: Option<Vec<u32>>, // Σ 1 over non-blacklisted positions

    cov_stage: Stage,  // Lifecycle for coverage
    bl_stage: BlStage, // Lifecycle for blacklist
}

impl CoveragePrefix {
    /// initialize_coverage_prefix
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
    pub fn initialize_coverage_prefix(length: u32) -> Self {
        Self {
            length,
            delta: vec![0.0; length as usize + 1],
            bl_delta: None,
            coverage: None,
            bl_mask: None,
            psum_all: None,
            psum_allowed: None,
            psum_allowed_count: None,
            cov_stage: Stage::Building,
            bl_stage: BlStage::Absent,
        }
    }

    /// add_fragment_to_prefix: +1 at start, -1 at end
    ///
    /// Parameters
    /// ----------
    /// - frag:
    ///     Fragment on the reference `[start, end)`, 0-based, end-exclusive.
    #[inline]
    pub fn add_fragment_to_prefix(&mut self, frag: Fragment) -> Result<()> {
        self.add_fragment_to_prefix_weighted(frag, 1.0)
    }

    /// add_fragment_to_prefix with floating weight w: +w at start, -w at end
    ///
    /// Parameters
    /// ----------
    /// - frag:
    ///     Fragment on the reference `[start, end)`.
    /// - weight:
    ///     Weight to add, must be finite and >= 0.
    #[inline]
    pub fn add_fragment_to_prefix_weighted(&mut self, frag: Fragment, weight: f32) -> Result<()> {
        if !self.prefix_available() {
            anyhow::bail!(
                "prefix was dropped; cannot add fragments. Rebuild or create a new CoveragePrefix"
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
        let start = frag.start as usize;
        let end = frag.end as usize;

        if start >= end {
            anyhow::bail!("fragment start {} >= end {}", frag.start, frag.end);
        }
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

    /// finalize_coverage: build per-base coverage from the +w/-w prefix.
    ///
    /// Returns
    /// -------
    /// - coverage:
    ///     Borrowed slice of per-base coverage with length = `length`.
    pub fn finalize_coverage(&mut self) -> &[f32] {
        // Build coverage from the +w/-w prefix without destroying delta
        let mut cov = Vec::with_capacity(self.length as usize);
        cov.resize(self.length as usize, 0.0);

        // Cumulative sum over delta
        let mut run = 0.0_f64;
        for i in 0..=self.length as usize {
            run += self.delta[i] as f64; // f64 to reduce rounding error accumulation
            if i < self.length as usize {
                cov[i] = run as f32;
            }
        }

        self.coverage = Some(cov);
        self.invalidate_indexes();
        self.cov_stage = Stage::Covered;

        self.coverage.as_ref().unwrap()
    }

    /// initialize_blacklist_prefix: enable optional blacklist delta (+1/-1)
    pub fn initialize_blacklist_prefix(&mut self) {
        if self.bl_delta.is_none() {
            self.bl_delta = Some(vec![0; self.length as usize + 1]);
        }
        self.bl_stage = BlStage::Building;
    }

    /// add_blacklist_to_prefix: add a blacklist interval `[start, end)`
    ///
    /// Parameters
    /// ----------
    /// - start:
    ///     Start of interval, inclusive.
    /// - end:
    ///     End of interval, exclusive.
    pub fn add_blacklist_to_prefix(&mut self, start: u64, end: u64) -> Result<()> {
        self.initialize_blacklist_prefix();

        {
            // Limit the mutable borrow of bl_delta to this block
            let bl = self
                .bl_delta
                .as_mut()
                .expect("blacklist delta must be initialized");
            let n = bl.len();
            let s = start as usize;
            let e = end as usize;

            if s >= e {
                anyhow::bail!("blacklist start {} >= end {}", start, end);
            }
            if e > self.length as usize || e >= n {
                anyhow::bail!(
                    "blacklist end {} out of bounds for sequence length {}",
                    end,
                    self.length
                );
            }

            bl[s] = bl[s].saturating_add(1);
            bl[e] = bl[e].saturating_sub(1);
        } // bl borrow ends here

        // Editing the prefix invalidates any finalized mask and indexes
        self.bl_mask = None;
        self.bl_stage = BlStage::Building;
        self.invalidate_indexes();

        Ok(())
    }

    /// add_blacklist_many_to_prefix: add multiple blacklist intervals `[start, end)` in one pass
    ///
    /// Parameters
    /// ----------
    /// - intervals:
    ///     Intervals as pairs in `[start, end)`; `start < end` and `end <= length`
    pub fn add_blacklist_many_to_prefix(&mut self, intervals: &[(u64, u64)]) -> Result<()> {
        self.initialize_blacklist_prefix();

        {
            // Limit mutable borrow of bl_delta to this block
            let bl = self
                .bl_delta
                .as_mut()
                .expect("blacklist delta must be initialized");

            let n = bl.len();
            let len_u64 = self.length as u64;

            for &(s64, e64) in intervals {
                if s64 >= e64 {
                    anyhow::bail!("blacklist start {} >= end {}", s64, e64);
                }
                if e64 > len_u64 {
                    anyhow::bail!(
                        "blacklist end {} out of bounds for sequence length {}",
                        e64,
                        self.length
                    );
                }

                let s = s64 as usize;
                let e = e64 as usize;

                if e >= n {
                    anyhow::bail!(
                        "blacklist end {} out of bounds for internal prefix length {}",
                        e,
                        n
                    );
                }

                bl[s] = bl[s].saturating_add(1);
                bl[e] = bl[e].saturating_sub(1);
            }
        }

        // Edits invalidate any finalized mask and indexes
        self.bl_mask = None;
        self.bl_stage = BlStage::Building;
        self.invalidate_indexes();
        Ok(())
    }

    /// finalize_blacklist_prefix: convert +1/-1 delta -> per-base mask where 1 = blacklisted, 0 = allowed
    pub fn finalize_blacklist_prefix(&mut self) {
        let Some(bl) = self.bl_delta.as_ref() else {
            // CHANGED: as_ref()
            self.bl_mask = None;
            self.invalidate_indexes();
            self.bl_stage = BlStage::Absent; // NEW
            return;
        };

        // Cumulative sum over blacklist delta
        let mut run = 0i32;
        let mut mask = vec![0u8; self.length as usize];
        // Walk bl_delta non-destructively
        for i in 0..=self.length as usize {
            run += bl[i];
            if i < mask.len() {
                mask[i] = if run > 0 { 1 } else { 0 };
            }
        }
        self.bl_mask = Some(mask);
        // Invalidate indexes to rebuild masked aggregates
        self.invalidate_indexes();
        self.bl_stage = BlStage::Finalized;
    }

    /// build_query_index: prepare prefix sums for fast interval queries
    ///
    /// Returns
    /// -------
    /// - _:
    ///     Err if coverage has not been finalized.
    pub fn build_query_index(&mut self) -> Result<()> {
        // Ensure per-base coverage is available
        let cov = match self.coverage.as_ref() {
            Some(c) => c,
            None => {
                anyhow::bail!("coverage not finalized, call finalize_coverage() first")
            }
        };

        // Number of positions
        let n = cov.len();

        // Allocate prefix arrays of length n+1
        // Index 0 stores the empty prefix so sums over [a, b) are psum[b] - psum[a]
        let mut psum_all = Vec::with_capacity(n + 1);
        let mut psum_allowed = Vec::with_capacity(n + 1);
        let mut psum_allowed_count = Vec::with_capacity(n + 1);

        // Prefix base case at index 0
        psum_all.push(0.0_f64);
        psum_allowed.push(0.0_f64);
        psum_allowed_count.push(0u32);

        match self.bl_mask.as_ref() {
            Some(mask) => {
                // Mask present -> build two parallel prefix sums
                // psum_all accumulates coverage at every base
                // psum_allowed accumulates coverage only at unmasked bases
                // psum_allowed_count counts unmasked bases for use as the average denominator
                for i in 0..n {
                    let c = cov[i] as f64;
                    let allowed = mask[i] == 0;

                    // Read previous prefix values
                    let prev_all = *psum_all.last().unwrap();
                    let prev_allow = *psum_allowed.last().unwrap();
                    let prev_cnt = *psum_allowed_count.last().unwrap();

                    // Always include coverage in psum_all
                    psum_all.push(prev_all + c);

                    // Include coverage and count only if allowed
                    psum_allowed.push(prev_allow + if allowed { c } else { 0.0 });
                    psum_allowed_count.push(prev_cnt + if allowed { 1 } else { 0 });
                }
            }
            None => {
                // No mask -> allowed sums equal all sums and the count is simply i+1
                for i in 0..n {
                    let c = cov[i] as f64;
                    let prev_all = *psum_all.last().unwrap();

                    // Update all-bases prefix sum
                    psum_all.push(prev_all + c);

                    // Reuse the newest psum_all value for psum_allowed
                    psum_allowed.push(*psum_all.last().unwrap());

                    // Every base is allowed so count increases by one
                    psum_allowed_count.push((i as u32) + 1);
                }
            }
        }

        // Store results on the struct
        self.psum_all = Some(psum_all);
        self.psum_allowed = Some(psum_allowed);
        self.psum_allowed_count = Some(psum_allowed_count);

        // Mark coverage stage as Indexed
        self.cov_stage = Stage::Indexed;

        Ok(())
    }

    /// sum_coverage: sum of coverage in `[start, end)`
    ///
    /// Parameters
    /// ----------
    /// - start:
    ///     Start of interval, inclusive.
    /// - end:
    ///     End of interval, exclusive.
    /// - exclude_blacklisted:
    ///     Exclude blacklisted positions from the sum.
    ///
    /// Returns
    /// -------
    /// - sum:
    ///     Coverage sum over the interval, masked if requested.
    pub fn sum_coverage(&mut self, start: u32, end: u32, exclude_blacklisted: bool) -> Result<f64> {
        self.ensure_coverage()?;
        self.ensure_indexes()?;
        self.ensure_mask_if_excluding(exclude_blacklisted)?;
        self.check_bounds(start, end)?;

        let a = start as usize;
        let b = end as usize;

        let sum = if exclude_blacklisted {
            let pa = self.psum_allowed.as_ref().unwrap();
            pa[b] - pa[a]
        } else {
            let pa = self.psum_all.as_ref().unwrap();
            pa[b] - pa[a]
        };
        Ok(sum)
    }

    /// avg_coverage: average coverage in `[start, end)`
    ///
    /// Parameters
    /// ----------
    /// - start:
    ///     Start of interval, inclusive.
    /// - end:
    ///     End of interval, exclusive.
    /// - exclude_blacklisted:
    ///     Exclude blacklisted positions from the average.
    ///
    /// Returns
    /// -------
    /// - avg:
    ///     Coverage average over the interval, masked if requested.
    pub fn avg_coverage(&mut self, start: u32, end: u32, exclude_blacklisted: bool) -> Result<f32> {
        self.ensure_coverage()?;
        self.ensure_indexes()?;
        self.ensure_mask_if_excluding(exclude_blacklisted)?;
        self.check_bounds(start, end)?;

        let a = start as usize;
        let b = end as usize;

        if exclude_blacklisted {
            let pa = self.psum_allowed.as_ref().unwrap();
            let cnt = self.psum_allowed_count.as_ref().unwrap();
            let sum = pa[b] - pa[a];
            let n_ok = (cnt[b] - cnt[a]) as u32;
            if n_ok == 0 {
                return Ok(0.0);
            }
            Ok((sum / n_ok as f64) as f32)
        } else {
            let pa = self.psum_all.as_ref().unwrap();
            let span = end - start;
            if span == 0 {
                return Ok(0.0);
            }
            let sum = pa[b] - pa[a];
            Ok((sum / span as f64) as f32)
        }
    }

    /// Batch sum over intervals using prefix sums.
    ///
    /// Parameters
    /// ----------
    /// - intervals:
    ///     Half-open intervals `[start, end)`.
    /// - exclude_blacklisted:
    ///     Exclude blacklisted positions from the sum.
    /// - parallelize:
    ///     Process intervals with rayon parallel iterators.
    ///
    /// Returns
    /// -------
    /// - sums:
    ///     A coverage sum per interval.
    pub fn bulk_sum_coverage(
        &mut self,
        intervals: &[(u32, u32)],
        exclude_blacklisted: bool,
        parallelize: bool,
    ) -> Result<Vec<f64>> {
        self.ensure_coverage()?;
        self.ensure_indexes()?;
        self.ensure_mask_if_excluding(exclude_blacklisted)?;
        for &(a, b) in intervals {
            self.check_bounds(a, b)?;
        }

        if exclude_blacklisted {
            let pa = self.psum_allowed.as_ref().unwrap();
            let f = |&(a, b): &(u32, u32)| -> f64 {
                let a = a as usize;
                let b = b as usize;
                pa[b] - pa[a]
            };
            Ok(if parallelize {
                intervals.par_iter().map(f).collect()
            } else {
                intervals.iter().map(f).collect()
            })
        } else {
            let pa = self.psum_all.as_ref().unwrap();
            let f = |&(a, b): &(u32, u32)| -> f64 {
                let a = a as usize;
                let b = b as usize;
                pa[b] - pa[a]
            };
            Ok(if parallelize {
                intervals.par_iter().map(f).collect()
            } else {
                intervals.iter().map(f).collect()
            })
        }
    }

    /// Batch average over intervals using prefix sums.
    ///
    /// Parameters
    /// ----------
    /// - intervals:
    ///     Half-open intervals `[start, end)`.
    /// - exclude_blacklisted:
    ///     Exclude blacklisted positions from the averages.
    /// - parallelize:
    ///     Process intervals with rayon parallel iterators.
    ///
    /// Returns
    /// -------
    /// - avgs:
    ///     An average coverage per interval.
    pub fn bulk_avg_coverage(
        &mut self,
        intervals: &[(u32, u32)],
        exclude_blacklisted: bool,
        parallelize: bool,
    ) -> Result<Vec<f32>> {
        self.ensure_coverage()?;
        self.ensure_indexes()?;
        self.ensure_mask_if_excluding(exclude_blacklisted)?;
        for &(a, b) in intervals {
            self.check_bounds(a, b)?;
        }

        if exclude_blacklisted {
            let pa = self.psum_allowed.as_ref().unwrap();
            let cnt = self.psum_allowed_count.as_ref().unwrap();
            let f = |&(a, b): &(u32, u32)| -> f32 {
                let a = a as usize;
                let b = b as usize;
                let sum = pa[b] - pa[a];
                let n_ok = (cnt[b] - cnt[a]) as u32;
                if n_ok == 0 {
                    0.0
                } else {
                    (sum / n_ok as f64) as f32
                }
            };
            Ok(if parallelize {
                intervals.par_iter().map(f).collect()
            } else {
                intervals.iter().map(f).collect()
            })
        } else {
            let pa = self.psum_all.as_ref().unwrap();
            let f = |&(a, b): &(u32, u32)| -> f32 {
                let a = a as usize;
                let b = b as usize;
                let span = b - a;
                if span == 0 {
                    0.0
                } else {
                    ((pa[b] - pa[a]) / span as f64) as f32
                }
            };
            Ok(if parallelize {
                intervals.par_iter().map(f).collect()
            } else {
                intervals.iter().map(f).collect()
            })
        }
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
    /// Value is 1 for blacklisted, 0 for allowed.
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

    /// coverage: borrowed per-base coverage slice if finalized
    ///
    /// Returns
    /// -------
    /// - coverage:
    ///     Per-base coverage of length `length` if available.
    pub fn coverage(&self) -> Option<&[f32]> {
        self.coverage.as_deref()
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

    /// length accessor
    pub fn length(&self) -> u32 {
        self.length
    }

    /// Drop the +w/-w delta to free memory. Further add_* calls will error.
    pub fn drop_prefix(&mut self) {
        self.delta.clear();
        self.delta.shrink_to_fit();
    }

    // Helper for testing
    pub fn _get_bl_delta(&self) -> Option<Vec<i32>> {
        self.bl_delta.clone()
    }

    // Private helpers

    fn ensure_indexes(&mut self) -> Result<()> {
        if self.psum_all.is_none()
            || self.psum_allowed.is_none()
            || self.psum_allowed_count.is_none()
        {
            self.build_query_index()?;
        }
        Ok(())
    }

    fn invalidate_indexes(&mut self) {
        self.psum_all = None;
        self.psum_allowed = None;
        self.psum_allowed_count = None;
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

    #[inline]
    fn prefix_available(&self) -> bool {
        !self.delta.is_empty()
    }

    #[inline]
    fn ensure_mask_if_excluding(&self, exclude_blacklisted: bool) -> Result<()> {
        if exclude_blacklisted {
            if matches!(self.bl_stage, BlStage::Building) {
                anyhow::bail!(
                    "blacklist present but not finalized; call finalize_blacklist_prefix()"
                );
            }
            // If `Absent`, we allow excluding (no mask means nothing to exclude).
        }
        Ok(())
    }
}
