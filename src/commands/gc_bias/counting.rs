use anyhow::{Result, ensure};
use ndarray::{Array2, ArrayBase, Axis, Data, DataMut, Ix2};

use crate::commands::gc_bias::smoothing::smooth_length_row_in_place;
use crate::shared::interval::{IndexedInterval, Interval};

/// Prefix sums (cumsum) to compute GC and AT fractions while excluding Ns.
/// `gc[i]`   = # of G/C in seq[0..i)
/// `acgt[i]`= # of A/T/G/C in seq[0..i)
pub struct GCPrefixes {
    pub gc: Vec<u32>,
    pub acgt: Vec<u32>,
}

impl GCPrefixes {
    /// Return the GC count in a checked half-open interval.
    #[inline]
    pub fn gc_count(&self, interval: Interval<usize>) -> Result<u32> {
        ensure!(
            interval.end() < self.gc.len(),
            "GC interval [{}, {}) out of bounds for prefix length {}",
            interval.start(),
            interval.end(),
            self.gc.len().saturating_sub(1)
        );
        Ok(self.gc[interval.end()] - self.gc[interval.start()])
    }

    /// Return the A/C/G/T count in a checked half-open interval.
    #[inline]
    pub fn acgt_count(&self, interval: Interval<usize>) -> Result<u32> {
        ensure!(
            interval.end() < self.acgt.len(),
            "ACGT interval [{}, {}) out of bounds for prefix length {}",
            interval.start(),
            interval.end(),
            self.acgt.len().saturating_sub(1)
        );
        Ok(self.acgt[interval.end()] - self.acgt[interval.start()])
    }
}

/// Build prefix-sums (cumsum) for GC and ACGT (non-N) counts over a byte slice.
/// This lets you compute GC% on A/T/G/C only, so AT% = 1 - GC%.
///
/// Ignores Ns.
pub fn build_gc_prefixes(seq: &[u8]) -> GCPrefixes {
    let mut gc = Vec::with_capacity(seq.len() + 1);
    let mut acgt = Vec::with_capacity(seq.len() + 1);
    gc.push(0);
    acgt.push(0);

    for &b in seq {
        let is_gc = matches!(b, b'G' | b'g' | b'C' | b'c') as u32;
        let is_acgt = matches!(b, b'A' | b'a' | b'T' | b't' | b'G' | b'g' | b'C' | b'c') as u32;

        gc.push(gc.last().copied().unwrap() + is_gc);
        acgt.push(acgt.last().copied().unwrap() + is_acgt);
    }

    GCPrefixes { gc, acgt }
}

/// Compute the GC integer percentage for a window [start, end), excluding 'N's.
///
/// `min_acgt_count`: Minimum number of actual ACGT bases counted in the window.
///   E.g. if most of the window is blacklisted or Ns.
///
/// Returns `None` if the window has too few A/T/G/C bases.
#[inline]
pub fn get_gc_integer_percentage_for_window(
    prefixes: &GCPrefixes,
    window: Interval<usize>,
    min_acgt_fraction: f32,
    min_acgt_count: u32,
) -> Result<Option<usize>> {
    let gc = prefixes.gc_count(window)?;
    let acgt = prefixes.acgt_count(window)?;

    if acgt == 0 || acgt < min_acgt_count || (acgt as f32 / window.len() as f32) < min_acgt_fraction
    {
        return Ok(None);
    }

    // Use the same integer rounding as the ref-gc-bias tool!
    let gc_percent_bin = calculate_gc_bin(gc as u64, acgt as u64);
    Ok(Some(gc_percent_bin))
}

/// Count matrix for fragment coverage across GC fraction bins and fragment lengths.
///
/// While counting, the counts matrix is flattened (for contiguous memory).
/// Use `.as_array2()` to get a 2-dimensional Array where:
/// - Rows correspond to fragment lengths .
/// - Columns correspond to GC fraction bins.
#[derive(Debug, Clone)]
pub struct GCCounts {
    counts: Vec<f64>,
    pub length_min: usize,
    pub length_max: usize,
    num_lengths: usize,
    offsets: Vec<usize>,
    end_offset: usize,
    pub num_acgt_out_of: (u64, u64),
}

impl GCCounts {
    /// Create a new `GCCounts` with specified ranges and binning.
    ///
    /// Parameters
    /// ----------
    /// length_min: usize
    ///     Minimum fragment length (inclusive).
    /// length_max: usize
    ///     Maximum fragment length (inclusive).
    /// end_offset: usize
    ///     Number of bases trimmed from each fragment end when counting GC.
    /// num_acgt_out_of:
    ///     Number of ACGT bases and the total number of positions.
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     A `GCCounts` object with all counts initialized to zero.
    pub fn new(
        length_min: usize,
        length_max: usize,
        end_offset: usize,
        num_acgt_out_of: (u64, u64),
    ) -> Result<Self> {
        let (counts, offsets, num_lengths) =
            Self::initialize_counts(length_min, length_max, end_offset)?;
        Ok(Self {
            counts,
            length_min,
            length_max,
            num_lengths,
            offsets,
            end_offset,
            num_acgt_out_of,
        })
    }

    /// Initialize the counts with 0s and create the offsets for length-indexing.
    fn initialize_counts(
        length_min: usize,
        length_max: usize,
        end_offset: usize,
    ) -> Result<(Vec<f64>, Vec<usize>, usize)> {
        ensure!(length_min <= length_max, "length_min must be <= length_max");
        let num_lengths = length_max - length_min + 1;
        let mut offsets = Vec::with_capacity(num_lengths + 1);
        let mut acc: usize = 0;
        for length in length_min..=length_max {
            offsets.push(acc);
            let effective_length = length.saturating_sub(2 * end_offset);
            acc += effective_length + 1; // gc = 0..=effective_length
        }
        offsets.push(acc); // End index
        let counts = vec![0.0; acc];
        Ok((counts, offsets, num_lengths))
    }

    /// Check whether `(length, gc)` is within configured ranges.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length to test.
    /// gc: usize
    ///     GC count to test.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     True if both indices are in range.
    #[inline]
    fn in_bounds(&self, length: usize, gc: usize) -> bool {
        self.effective_length(length)
            .is_some_and(|effective_length| gc <= effective_length)
    }

    #[inline]
    fn effective_length(&self, length: usize) -> Option<usize> {
        if self.length_range().contains(&length) {
            Some(length.saturating_sub(2 * self.end_offset))
        } else {
            None
        }
    }

    /// Get index in the raw (flattened) counts for a given fragment length and GC count.
    #[inline]
    fn flat_index(&self, length: usize, gc: usize) -> Option<usize> {
        if self.in_bounds(length, gc) {
            let row = length - self.length_min;
            let start = self.offsets[row];
            Some(start + gc)
        } else {
            None
        }
    }

    /// Increment the counter for a given fragment length and GC count.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC count (number of GC bases after end offsets), absolute.
    pub fn incr(&mut self, length: usize, gc: usize) {
        if let Some(idx) = self.flat_index(length, gc) {
            self.counts[idx] += 1.0;
        }
    }

    /// Increment the counter by a weight for a given fragment length and GC count.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC count (number of GC bases after end offsets), absolute.
    /// weight: f64
    ///     Weight to count up.
    pub fn incr_weighted(&mut self, length: usize, gc: usize, weight: f64) {
        if let Some(idx) = self.flat_index(length, gc) {
            self.counts[idx] += weight;
        }
    }

    /// Get the count at a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    ///
    /// Returns
    /// -------
    /// count: Option<f64>
    ///     The count if indices are in range, otherwise `None`.
    #[allow(dead_code)]
    pub fn get(&self, length: usize, gc: usize) -> Option<f64> {
        self.flat_index(length, gc).map(|idx| self.counts[idx])
    }

    /// Set the count at a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    /// count: f64
    ///     Value to set as count.
    #[allow(dead_code)]
    pub fn set(&mut self, length: usize, gc: usize, count: f64) {
        if let Some(idx) = self.flat_index(length, gc) {
            self.counts[idx] = count;
        }
    }

    /// Borrow the raw counts (flattened vector). Non-mutable borrow.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    ///
    /// Returns
    /// -------
    /// count: Option<f64>
    ///     The count if indices are in range, otherwise `None`.
    pub fn borrow_raw_counts(&self) -> &Vec<f64> {
        &self.counts
    }

    /// Number of length rows.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of length bins.
    #[inline]
    pub fn n_lengths(&self) -> usize {
        self.num_lengths
    }

    /// Get the percentage of `ACGT` bases in the observed positions.
    /// I.e., the percentage of positions that are not blacklisted or ambiguous (`N`).
    ///
    /// Returns
    /// -------
    /// percentage: f64
    ///     Number of observed ACGT bases divided by total number of observed bases.
    ///     When no positions are observed, it returns `f64::NAN`.
    pub fn pct_acgt(&self) -> f64 {
        if self.num_acgt_out_of.1 == 0 {
            f64::NAN
        } else {
            100.0 * (self.num_acgt_out_of.0 as f64 / self.num_acgt_out_of.1 as f64)
        }
    }

    /// Sum of all counts.
    pub fn sum(&self) -> f64 {
        self.counts.iter().copied().sum()
    }

    /// Mean count across all GC-by-length cells. Returns 0.0 if empty.
    pub fn mean(&self) -> f64 {
        if self.counts.is_empty() {
            0.0
        } else {
            self.sum() / self.counts.len() as f64
        }
    }

    /// Get the flat slice bounds for a specific fragment length.
    pub fn length_bounds(&self, length: usize) -> Result<(usize, usize)> {
        ensure!(
            self.length_range().contains(&length),
            "length {} outside [{}, {}]",
            length,
            self.length_min,
            self.length_max
        );
        let row_idx = length - self.length_min;
        Ok((self.offsets[row_idx], self.offsets[row_idx + 1]))
    }

    /// Sum counts for a specific fragment length row.
    pub fn sum_for_length(&self, length: usize) -> Result<f64> {
        let (start, end) = self.length_bounds(length)?;
        Ok(self.counts[start..end].iter().copied().sum())
    }

    /// Scale all counts by `factor` in place.
    pub fn scale_counts(&mut self, factor: f64) -> Result<()> {
        ensure!(factor > 0.0);
        if factor == 1.0 {
            return Ok(());
        }
        for v in self.counts.iter_mut() {
            *v *= factor;
        }
        Ok(())
    }

    /// Create a zero-initialized `GCCounts` with the same ranges and shape as `self`.
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     New object with all counts set to zero and identical configuration.
    #[inline]
    pub fn zeroed_like(&self) -> Result<Self> {
        Self::new(self.length_min, self.length_max, self.end_offset, (0, 0))
    }

    /// Clear counts in place without reallocating.
    #[inline]
    pub fn clear(&mut self) {
        self.counts.fill(0.0);
        self.num_acgt_out_of = (0, 0);
    }

    #[inline]
    pub fn buffer_len(&self) -> usize {
        self.counts.len()
    }

    /// Build a `GCCounts` from raw components.
    ///
    /// Parameters
    /// ----------
    /// counts: Vec<f64>
    ///     Flattened counts buffer matching the internal layout for the provided ranges.
    /// length_min: usize
    ///     Minimum fragment length (inclusive).
    /// length_max: usize
    ///     Maximum fragment length (inclusive).
    /// end_offset: usize
    ///     Number of bases trimmed from each fragment end when counting GC.
    /// num_acgt_out_of:
    ///     Number of ACGT bases and the total number of positions.
    pub fn from_parts(
        counts: Vec<f64>,
        length_min: usize,
        length_max: usize,
        end_offset: usize,
        num_acgt_out_of: (u64, u64),
    ) -> Result<Self> {
        let (_init_counts, offsets, num_lengths) =
            Self::initialize_counts(length_min, length_max, end_offset)?;
        ensure!(
            counts.len() == offsets.last().copied().unwrap_or(0),
            "counts buffer length {} does not match expected {}",
            counts.len(),
            offsets.last().copied().unwrap_or(0)
        );
        Ok(Self {
            counts,
            length_min,
            length_max,
            num_lengths,
            offsets,
            end_offset,
            num_acgt_out_of,
        })
    }

    /// Check if two `GCCounts` are compatible for merging (same ranges and shape).
    ///
    /// Parameters
    /// ----------
    /// other: &GCCounts
    ///     The other counts matrix to compare with.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     `true` if the two objects have identical ranges and matrix shape.
    #[inline]
    pub fn is_compatible_with(&self, other: &Self) -> bool {
        self.length_min == other.length_min
            && self.length_max == other.length_max
            && self.n_lengths() == other.n_lengths()
            && self.offsets == other.offsets
            && self.counts.len() == other.counts.len()
            && self.end_offset == other.end_offset
    }

    /// Merge (sum) counts from `other` into `self` using saturating addition.
    ///
    /// Parameters
    /// ----------
    /// other: &GCCounts
    ///     Counts to add into `self`. Must be compatible (same ranges/shape).
    ///
    /// Returns
    /// -------
    /// result: Result<(), anyhow::Error>
    ///     An error is returned if the two objects are incompatible.
    pub fn merge_from(&mut self, other: &Self) -> Result<()> {
        if !self.is_compatible_with(other) {
            return Err(anyhow::anyhow!(
                "incompatible GCCounts: self={} vs other={}",
                self,
                other
            ));
        }
        for (idx, other_count) in other.borrow_raw_counts().iter().enumerate() {
            self.counts[idx] += other_count;
        }
        self.num_acgt_out_of.0 += other.num_acgt_out_of.0;
        self.num_acgt_out_of.1 += other.num_acgt_out_of.1;

        Ok(())
    }

    /// Collapse (sum) an iterator of `GCCounts` into a single object.
    ///
    /// All inputs must be compatible (same ranges/shape). Uses saturating addition.
    ///
    /// Parameters
    /// ----------
    /// iter: IntoIterator<Item = &GCCounts>
    ///     Collection of references to `GCCounts` to be summed.
    ///
    /// Returns
    /// -------
    /// total: GCCounts
    ///     The element-wise sum across all inputs.
    ///
    /// Examples
    /// --------
    /// ```ignore
    /// // Sum across chromosomes:
    /// // let by_chr: HashMap<String, GCCounts> = ...;
    /// let total = GCCounts::collapse(by_chr.values())?;
    /// ```
    #[allow(dead_code)]
    pub fn collapse<'a, I>(iter: I) -> Result<Self>
    where
        I: IntoIterator<Item = &'a GCCounts>,
    {
        let mut it = iter.into_iter();
        let first = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("collapse requires at least one GCCounts"))?;
        let mut acc = first.clone(); // Start from the first
        for g in it {
            acc.merge_from(g)?; // Then merge the rest
        }
        Ok(acc)
    }

    /// Slice the row for `length` from a flat counts buffer and smooth it in place.
    pub fn smooth_length_rows_in_place(&mut self, sigma: f64, radius: u8) -> Result<()> {
        for length in self.length_range() {
            let Some(effective_length) = self.effective_length(length) else {
                continue;
            };
            if effective_length == 0 {
                continue;
            }
            // Uses offsets based on original lengths
            // so end_offsets are already considered
            smooth_length_row_in_place(
                &mut self.counts,
                &self.offsets,
                self.length_min,
                length,
                sigma,
                radius as usize,
            )?;
        }
        Ok(())
    }

    pub fn length_range(&self) -> std::ops::RangeInclusive<usize> {
        self.length_min..=self.length_max
    }

    #[inline]
    pub fn end_offset(&self) -> usize {
        self.end_offset
    }

    /// Collapse ragged GC-count rows into a rectangular GC% x length matrix.
    ///
    /// The stored `end_offset` is applied to derive effective lengths before
    /// computing GC percentages, so bins that are theoretically unobservable
    /// are never populated.
    pub fn to_gc_percent_grid(&self, gc_min: usize, gc_max: usize) -> Result<Array2<f64>> {
        ensure!(gc_min < gc_max, "gc_min must be < gc_max");
        let num_lengths = self.length_max - self.length_min + 1;
        let num_gc_bins = gc_max - gc_min + 1;
        let mut grid = Array2::<f64>::zeros((num_lengths, num_gc_bins));
        for length in self.length_range() {
            let row_idx = length - self.length_min;
            let (start, end) = self.length_bounds(length)?;
            let effective_length = length.saturating_sub(2 * self.end_offset);
            if effective_length == 0 {
                continue;
            }
            for (gc, &val) in (0..).zip(self.counts[start..end].iter()) {
                let gc_pct = calculate_gc_bin(gc as u64, effective_length as u64);
                if (gc_min..=gc_max).contains(&gc_pct) {
                    let col_idx = gc_pct - gc_min;
                    grid[(row_idx, col_idx)] += val;
                }
            }
        }
        Ok(grid)
    }
}

impl std::fmt::Display for GCCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GCCounts(len:[{}..={}], pct_acgt:({}) )",
            self.length_min,
            self.length_max,
            self.pct_acgt(),
        )
    }
}

/// Count reference GC per fragment length for every window on one chromosome.
///
/// For each window `[start, end)` and each given start position `s` within that
/// window (`start <= s < end`), this function considers all fragment lengths
/// `L` in `[min_len, max_len)` such that `s + L <= end`. It uses prefix arrays to
/// compute, in O(1), both:
/// - the number of A/C/G/T bases (excluding Ns/blacklist) in `[s, s+L)`, and
/// - the number of G or C bases in `[s, s+L)`.
///
/// A window/length is **counted** only if it meets both ACGT requirements:
/// - `acgt_count >= min_acgt_count`, and
/// - `acgt_count / L >= min_acgt_fraction`.
///
/// When counted, the raw GC **count** for the fragment's effective length
/// (after applying `end_offset`) is stored. Percentages are derived later with
/// `to_gc_percent_grid`, which uses the effective length for the same
/// half-up rounding as `calculate_gc_bin`.
///
/// The sampled `start_positions` must be **sorted, unique**, and refer to the same
/// chromosome as `gc_prefixes`. Windows are assumed **sorted by start** (they may
/// overlap). The function advances a pointer through `start_positions` as windows
/// progress to avoid re-scanning earlier starts.
///
/// Parameters
/// ----------
/// counts_by_bin: &mut Vec<GCCounts>
///     Per-window accumulator. For window index `i`, `counts_by_bin[i]` is updated
///     by calling `incr(fragment_length, gc_count)`.
/// gc_prefixes: &GCPrefixes
///     Prefix arrays with one extra sentinel element: `gc[k]` and `acgt[k]` give
///     cumulative counts in `[0, k)`. Requires `gc.len() == acgt.len() >= chrom_len + 1`.
/// length_range: (u64, u64)
///     Half-open `[min_len, max_len)` in base pairs for fragment lengths.
/// windows: &[(u64, u64, u64)]
///     Start-sorted windows as `(start, end, original_idx)`. Each window is clamped
///     to `chrom_len`.
/// start_positions: &[usize]
///     Sorted, unique genomic start indices for this chromosome.
/// chrom_len: u64
///     Chromosome length; used to cap window ends.
/// min_acgt_fraction: f32
///     Minimum fraction of A/C/G/T within the fragment (after masking Ns/blacklist).
/// min_acgt_count: u32
///     Minimum absolute count of A/C/G/T within the fragment.
///
/// Returns
/// -------
/// None
///     Updates `counts_by_bin` in place.
pub fn count_reference_gc_and_length_by_window(
    counts_by_bin: &mut Vec<GCCounts>,
    gc_prefixes: &GCPrefixes,
    length_range: (u64, u64), // [min_len, max_len) in bp
    windows: &[IndexedInterval<u64>],
    start_positions: &[usize], // Sorted unique genomic starts for this chromosome
    chrom_len: u64,
    min_acgt_fraction: f32, // E.g., 0.8
    min_acgt_count: u32,
    end_offset: usize,
) -> Result<()> {
    ensure!(
        gc_prefixes.gc.len() == gc_prefixes.acgt.len(),
        "GC and ACGT prefix lengths differed: {} vs {}",
        gc_prefixes.gc.len(),
        gc_prefixes.acgt.len()
    );
    ensure!(
        gc_prefixes.gc.len() >= chrom_len as usize + 1,
        "prefix arrays should be chrom_len+1"
    );

    let min_len = length_range.0 as usize;
    let max_len = length_range.1 as usize; // exclusive

    // Precompute required ACGT counts per length: max(ceil(frac * len), min_count)
    let mut required_acgt_per_len = vec![0u32; max_len + 1];
    for len in min_len..max_len {
        let effective_length = len.saturating_sub(2 * end_offset);
        if effective_length == 0 {
            continue;
        }
        let req_by_frac = (min_acgt_fraction * (effective_length as f32)).ceil() as u32;
        required_acgt_per_len[len] = req_by_frac.max(min_acgt_count).max(1);
    }

    // Pointer into `start_positions`, advanced monotonically as windows progress
    let mut start_ptr = 0usize;

    for (win_idx, window) in windows.iter().enumerate() {
        let window_start = window.start();
        let window_end = window.end().min(chrom_len);
        let window_end_usize = window_end as usize;
        let window_start_usize = window_start as usize;
        let window_len = window_end_usize - window_start_usize;
        if window_len == 0 {
            continue;
        }

        // Advance the start pointer to the first start >= window_start.
        while start_ptr < start_positions.len() && start_positions[start_ptr] < window_start_usize {
            start_ptr += 1;
        }

        // Iterate all sampled starts that fall inside [window_start, window_end).
        let mut j = start_ptr;
        while j < start_positions.len() {
            let start_pos = start_positions[j];
            if start_pos >= window_end_usize {
                break; // Past the window
            }

            // Remaining room to the right edge of the window (inclusive end condition)
            // Valid fragment lengths satisfy: min_len <= frag_len <= (window_end - start_pos)
            let remaining = window_end_usize - start_pos;
            if remaining >= min_len {
                let frag_max_excl = max_len.min(remaining + 1);

                for frag_len in min_len..frag_max_excl {
                    // Apply end offsets to GC window
                    if frag_len <= 2 * end_offset {
                        continue;
                    }
                    let fragment_interval = Interval::new(
                        start_pos,
                        start_pos
                            .checked_add(frag_len)
                            .expect("fragment interval end should not overflow usize"),
                    )
                    .expect("fragment interval should stay non-empty");
                    let gc_interval = fragment_interval
                        .contract(end_offset)
                        .expect("GC interval should stay non-empty after end offsets");

                    // Prefix lookups
                    let acgt_count = gc_prefixes.acgt_count(gc_interval)?;
                    if acgt_count < required_acgt_per_len[frag_len] {
                        continue;
                    }

                    let gc_count = gc_prefixes.gc_count(gc_interval)?;

                    counts_by_bin[win_idx].incr(frag_len, gc_count as usize);
                }
            }

            j += 1;
        }
    }
    Ok(())
}

/// Round to the nearest percent (**half-up**) using integer math.
/// Integer division floors: (100*gc)/acgt would always round down (bias!).
/// Trick: add half the denominator before dividing -> values with fractional part >= 0.5 round up.
/// Formula: round(x/y) = (x + y/2) / y, for y > 0.
/// Here: x = 100 * gc_count, y = acgt_count.
/// Examples:
///  - gc=1, acgt=3 -> exact 33.33…% -> (100 + 1)/3 = 33
///  - gc=2, acgt=3 -> exact 66.66…% -> (200 + 1)/3 = 67
///  - gc=3, acgt=3 -> exact 100%     -> (300 + 1)/3 = 100 (then clamped to ≤100 below)
pub fn calculate_gc_bin(gc_count: u64, acgt_count: u64) -> usize {
    ((gc_count * 100 + (acgt_count / 2)) / acgt_count).min(100) as usize
}

/// Precompute how many raw GC counts map to each integer GC% bin for every length.
///
/// Used to debias the uneven bin widths caused by half-up rounding of `gc_count / length`
/// to an integer percent. Each output row corresponds to a fragment length and each
/// column is a GC% bin in `[0, 100]`. Entries store how many distinct `gc_count` values
/// round into that bin for the given length.
pub fn gc_percent_widths(length_min: usize, length_max: usize, end_offset: usize) -> Array2<u16> {
    assert!(
        length_max >= length_min,
        "length range must be non-empty ({}..={})",
        length_min,
        length_max
    );
    let num_lengths = length_max - length_min + 1;
    let mut widths = Array2::<u16>::zeros((num_lengths, 101));

    for length in length_min..=length_max {
        let row_idx = length - length_min;
        let effective_length = length.saturating_sub(2 * end_offset);
        if effective_length == 0 {
            continue;
        }
        for gc_count in 0..=effective_length {
            let gc_bin = calculate_gc_bin(gc_count as u64, effective_length as u64);
            // Widths are small (<= length+1), so u16 is sufficient.
            let current = widths
                .get_mut((row_idx, gc_bin))
                .expect("bin index in bounds");
            *current = current.saturating_add(1);
        }
    }

    widths
}

pub fn apply_gc_percent_width_correction<S, W>(
    counts: &mut ArrayBase<S, Ix2>,
    widths: &ArrayBase<W, Ix2>,
) -> Result<()>
where
    S: DataMut<Elem = f64>,
    W: Data<Elem = u16>,
{
    ensure!(
        counts.dim() == widths.dim(),
        "GC percent widths shape {:?} must match counts shape {:?}",
        widths.dim(),
        counts.dim()
    );

    let (n_rows, n_cols) = counts.dim();
    for row in 0..n_rows {
        let sum_before: f64 = counts.row(row).sum();
        let mut sum_after = 0.0;
        for col in 0..n_cols {
            let width = widths[(row, col)];
            if width > 0 {
                let val = counts[(row, col)] / width as f64;
                counts[(row, col)] = val;
                sum_after += val;
            } else {
                counts[(row, col)] = 0.0;
            }
        }
        if sum_after > 0.0 && sum_before > 0.0 {
            let scale = sum_before / sum_after;
            let mut row_view = counts.index_axis_mut(Axis(0), row);
            row_view *= scale;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    include!("counting_tests.rs");
}
