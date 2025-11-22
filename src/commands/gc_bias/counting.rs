use core::f64;

use ndarray::{Array2, Array3, s};

/// Prefix sums (cumsum) to compute GC and AT fractions while excluding Ns.
/// `gc[i]`   = # of G/C in seq[0..i)
/// `acgt[i]`= # of A/T/G/C in seq[0..i)
pub struct GCPrefixes {
    pub gc: Vec<u32>,
    pub acgt: Vec<u32>,
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
    start: usize,
    end: usize,
    min_acgt_fraction: f32,
    min_acgt_count: u32,
) -> Option<usize> {
    debug_assert!(
        start < end && end <= prefixes.gc.len() - 1,
        "GC window [{}, {}) out of bounds (len={})",
        start,
        end,
        prefixes.gc.len() - 1
    );

    let gc = prefixes.gc[end] - prefixes.gc[start];
    let acgt = prefixes.acgt[end] - prefixes.acgt[start];
    let length = end - start;

    if acgt == 0 || acgt < min_acgt_count || (acgt as f32 / length as f32) < min_acgt_fraction {
        return None;
    }

    // Use the same integer rounding as the reference-gc tool!
    let gc_percent_bin = calculate_gc_bin(gc as u64, acgt as u64);
    Some(gc_percent_bin)
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
    pub gc_min: usize,
    pub gc_max: usize,
    pub length_min: usize,
    pub length_max: usize,
    num_lengths: usize,
    num_gc_bins: usize,
    pub num_acgt_out_of: (u64, u64),
}

impl GCCounts {
    /// Create a new `GCCounts` with specified ranges and binning.
    ///
    /// Parameters
    /// ----------
    /// gc_min: usize
    ///     Minimum GC bin (inclusive).
    /// gc_max: usize
    ///     Maximum GC bin (inclusive).
    /// length_min: usize
    ///     Minimum fragment length (inclusive).
    /// length_max: usize
    ///     Maximum fragment length (inclusive).
    /// num_acgt_out_of:
    ///     Number of ACGT bases and the total number of positions.
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     A `GCCounts` object with all counts initialized to zero.
    pub fn new(
        gc_min: usize,
        gc_max: usize,
        length_min: usize,
        length_max: usize,
        num_acgt_out_of: (u64, u64),
    ) -> Self {
        let num_gc_bins = gc_max - gc_min + 1;
        let num_lengths = length_max - length_min + 1;
        let counts = vec![0f64; num_gc_bins * num_lengths];
        Self {
            counts,
            gc_min,
            gc_max,
            length_min,
            length_max,
            num_lengths,
            num_gc_bins,
            num_acgt_out_of,
        }
    }

    /// Check whether `(length, gc)` is within configured ranges.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length to test.
    /// gc: usize
    ///     GC bin to test.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     True if both indices are in range.
    #[inline]
    fn in_bounds(&self, length: usize, gc: usize) -> bool {
        (self.length_min..=self.length_max).contains(&length)
            && (self.gc_min..=self.gc_max).contains(&gc)
    }

    /// Compute row/column indices in `Array2` (see `.as_array2()`)
    /// output for `(length, gc)` if in bounds.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute, not zero-based).
    /// gc: usize
    ///     GC bin (absolute, not zero-based).
    ///
    /// Returns
    /// -------
    /// idx: Option<(usize, usize)>
    ///     `(row, col)` zero-based indices if in range, otherwise `None`.
    #[inline]
    pub fn index_of(&self, length: usize, gc: usize) -> Option<(usize, usize)> {
        if self.in_bounds(length, gc) {
            Some((length - self.length_min, gc - self.gc_min))
        } else {
            None
        }
    }

    /// Get index in the raw (flattened) counts for a given fragment length and gc.
    #[inline]
    pub fn flat_index(&self, length: usize, gc: usize) -> Option<usize> {
        if self.in_bounds(length, gc) {
            let row = length - self.length_min;
            let col = gc - self.gc_min;
            Some(row * self.num_gc_bins + col)
        } else {
            None
        }
    }

    /// Increment the counter for a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    pub fn incr(&mut self, length: usize, gc: usize) {
        if let Some(idx) = self.flat_index(length, gc) {
            self.counts[idx] += 1.0;
        }
    }

    /// Increment the counter by a weight for a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
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

    /// Get the counts as an `Array2` object.
    ///
    /// Returns
    /// -------
    /// counts: Array2<f64>
    pub fn to_array2(&self) -> Array2<f64> {
        Array2::from_shape_vec((self.num_lengths, self.num_gc_bins), self.counts.clone())
            .expect("inconsistent row lengths")
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

    /// Number of GC columns.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of GC bins.
    #[inline]
    pub fn n_gc_bins(&self) -> usize {
        self.num_gc_bins
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

    /// Create a zero-initialized `GCCounts` with the same ranges and shape as `self`.
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     New object with all counts set to zero and identical configuration.
    #[inline]
    pub fn zeroed_like(&self) -> Self {
        Self {
            counts: vec![0f64; self.num_gc_bins * self.num_lengths],
            gc_min: self.gc_min,
            gc_max: self.gc_max,
            length_min: self.length_min,
            length_max: self.length_max,
            num_gc_bins: self.num_gc_bins,
            num_lengths: self.num_lengths,
            num_acgt_out_of: (0, 0),
        }
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
        self.gc_min == other.gc_min
            && self.gc_max == other.gc_max
            && self.length_min == other.length_min
            && self.length_max == other.length_max
            && self.n_lengths() == other.n_lengths()
            && self.n_gc_bins() == other.n_gc_bins()
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
    pub fn merge_from(&mut self, other: &Self) -> anyhow::Result<()> {
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
    pub fn collapse<'a, I>(iter: I) -> anyhow::Result<Self>
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
}

impl Default for GCCounts {
    /// Create an empty default `GCCounts` (0–100 GC, 20–600 length).
    fn default() -> Self {
        Self::new(0, 100, 20, 600, (0, 0))
    }
}

impl std::fmt::Display for GCCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GCCounts(gc:[{}..={}], len:[{}..={}], dims:({},{}), pct_acgt:({}) )",
            self.gc_min,
            self.gc_max,
            self.length_min,
            self.length_max,
            self.n_lengths(),
            self.n_gc_bins(),
            self.pct_acgt(),
        )
    }
}

/// Stack counts from vector of `Array2` windows to a single 3d array.
pub fn stack_gc_counts(all_counts: &[Array2<f64>]) -> Array3<f64> {
    let n = all_counts.len();
    assert!(n > 0, "stack_gc_counts requires at least one window");

    let rows = all_counts[0].nrows();
    let cols = all_counts[0].ncols();

    let mut stacked = Array3::<f64>::zeros((n, rows, cols));
    for (idx, window) in all_counts.iter().enumerate() {
        assert_eq!(window.nrows(), rows, "mismatched length bins at {}", idx);
        assert_eq!(window.ncols(), cols, "mismatched GC bins at {}", idx);
        stacked.slice_mut(s![idx, .., ..]).assign(window);
    }

    stacked
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
/// When counted, the GC fraction is binned to a **percent** in `[0, 100]` using
/// **half-up rounding** *without floats*, via:
/// `round(100 * gc / acgt) = (100*gc + acgt/2) / acgt`.
/// This avoids the systematic low bias of floor‐division.
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
///     by calling `incr(fragment_length, gc_percent_bin)`.
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
    windows: &[(u64, u64, u64)],
    start_positions: &[usize], // sorted unique genomic starts for this chromosome
    chrom_len: u64,
    min_acgt_fraction: f32, // e.g., 0.8
    min_acgt_count: u32,
) {
    let gc_prefix = &gc_prefixes.gc; // prefix sums of GC counts
    let acgt_prefix = &gc_prefixes.acgt; // prefix sums of A/C/G/T (non-N/non-blacklist)
    debug_assert_eq!(gc_prefix.len(), acgt_prefix.len());
    debug_assert!(
        gc_prefix.len() >= chrom_len as usize + 1,
        "prefix arrays should be chrom_len+1"
    );

    let min_len = length_range.0 as usize;
    let max_len = length_range.1 as usize; // exclusive

    // Precompute required ACGT counts per length: max(ceil(frac * len), min_count)
    let mut required_acgt_per_len = vec![0u32; max_len + 1];
    for len in min_len..max_len {
        let req_by_frac = (min_acgt_fraction * (len as f32)).ceil() as u32;
        required_acgt_per_len[len] = req_by_frac.max(min_acgt_count).max(1);
    }

    // Pointer into `start_positions`, advanced monotonically as windows progress
    let mut start_ptr = 0usize;

    for (win_idx, &(window_start, mut window_end, _)) in windows.iter().enumerate() {
        window_end = window_end.min(chrom_len);
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
                let frag_max_excl = max_len.min(remaining as usize + 1);

                for frag_len in min_len..frag_max_excl {
                    let end_idx = start_pos + frag_len;

                    // Prefix lookups
                    let acgt_count = acgt_prefix[end_idx] - acgt_prefix[start_pos];
                    if acgt_count < required_acgt_per_len[frag_len] {
                        continue;
                    }

                    let gc_count = gc_prefix[end_idx] - gc_prefix[start_pos];

                    // Round to the nearest percent (**half-up**) using integer math.
                    let gc_percent_bin = calculate_gc_bin(gc_count as u64, acgt_count as u64);

                    counts_by_bin[win_idx].incr(frag_len, gc_percent_bin);
                }
            }

            j += 1;
        }
    }
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
    ((gc_count as u64 * 100 + (acgt_count as u64 / 2)) / acgt_count as u64).min(100) as usize
}
