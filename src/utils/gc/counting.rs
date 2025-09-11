use ndarray::Array3;

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

/// Compute the GC fraction for a window [start, end), excluding 'N's.
///
/// `min_acgt_count`: Minimum number of actual ACGT bases counted in the window.
///   E.g. if most of the window is blacklisted or Ns.
///
/// Returns None if the window has no A/T/G/C bases.
#[inline]
pub fn get_gc_fraction_in_window(
    prefixes: &GCPrefixes,
    start: usize,
    end: usize,
    min_acgt_fraction: f32,
    min_acgt_count: u32,
) -> Option<f32> {
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
    let gc_frac = gc as f32 / acgt as f32;
    Some(gc_frac)
}

/// Count matrix for fragment coverage across GC fraction bins and fragment lengths.
///
/// The matrix is two-dimensional:
/// - Rows correspond to fragment lengths.
/// - Columns correspond to GC fraction bins (0–100).
#[derive(Debug, Clone)]
pub struct GCCounts {
    pub counts: Vec<Vec<u64>>,
    pub gc_min: usize,
    pub gc_max: usize,
    pub length_min: usize,
    pub length_max: usize,
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
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     A `GCCounts` object with all counts initialized to zero.
    pub fn new(gc_min: usize, gc_max: usize, length_min: usize, length_max: usize) -> Self {
        let num_gc_bins = gc_max - gc_min + 1;
        let num_lengths = length_max - length_min + 1;
        let counts = vec![vec![0u64; num_gc_bins]; num_lengths];
        Self {
            counts,
            gc_min,
            gc_max,
            length_min,
            length_max,
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

    /// Compute row/column indices for `(length, gc)` if in bounds.
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

    /// Increment the counter for a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    pub fn incr(&mut self, length: usize, gc: usize) {
        if let Some((r, c)) = self.index_of(length, gc) {
            self.counts[r][c] = self.counts[r][c].saturating_add(1);
        }
    }

    // Get the count at a given fragment length and GC bin.
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
    /// count: Option<u32>
    ///     The count if indices are in range, otherwise `None`.
    pub fn get(&self, length: usize, gc: usize) -> Option<u64> {
        self.index_of(length, gc).map(|(r, c)| self.counts[r][c])
    }

    /// Number of length rows.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of length bins.
    #[inline]
    pub fn n_lengths(&self) -> usize {
        self.length_max - self.length_min + 1
    }

    /// Number of GC columns.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of GC bins.
    #[inline]
    pub fn n_gc_bins(&self) -> usize {
        self.gc_max - self.gc_min + 1
    }

    /// Create a zero-initialized `GCCounts` with the same ranges and shape as `self`.
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     New object with all counts set to zero and identical configuration.
    #[inline]
    pub fn zeroed_like(&self) -> Self {
        let n_len = self.n_lengths();
        let n_gc = self.n_gc_bins();
        Self {
            counts: vec![vec![0u64; n_gc]; n_len],
            gc_min: self.gc_min,
            gc_max: self.gc_max,
            length_min: self.length_min,
            length_max: self.length_max,
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
        for (r, row_other) in other.counts.iter().enumerate() {
            let row_self = &mut self.counts[r];
            for (c, &v) in row_other.iter().enumerate() {
                row_self[c] = row_self[c].saturating_add(v);
            }
        }
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
        Self::new(0, 100, 20, 600)
    }
}

impl std::fmt::Display for GCCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GCCounts(gc:[{}..={}], len:[{}..={}], dims:({},{}) )",
            self.gc_min,
            self.gc_max,
            self.length_min,
            self.length_max,
            self.n_lengths(),
            self.n_gc_bins()
        )
    }
}

/// Stack counts from vector of `GCCounts` to a single 3d array.
pub fn stack_gc_counts(all_counts: &Vec<GCCounts>) -> Array3<u64> {
    // Assume all GCCounts.counts have the same dimensions
    let n = all_counts.len();
    let l = all_counts[0].n_lengths();
    let g = all_counts[0].n_gc_bins();

    let mut arr = Array3::<u64>::zeros((n, l, g));
    for (i, gcc) in all_counts.iter().enumerate() {
        assert_eq!(gcc.n_lengths(), l, "mismatched length bins at {}", i);
        assert_eq!(gcc.n_gc_bins(), g, "mismatched GC bins at {}", i);
        for (li, row) in gcc.counts.iter().enumerate() {
            // row: Vec<u64> over GC bins
            for (gi, &val) in row.iter().enumerate() {
                arr[(i, li, gi)] = val;
            }
        }
    }

    arr
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
                    // Integer division floors: (100*gc)/acgt would always round down (bias!).
                    // Trick: add half the denominator before dividing -> values with fractional part ≥ 0.5 round up.
                    // Formula: round(x/y) = (x + y/2) / y, for y > 0.
                    // Here: x = 100 * gc_count, y = acgt_count.
                    // Examples:
                    //  - gc=1, acgt=3 -> exact 33.33…% -> (100 + 1)/3 = 33
                    //  - gc=2, acgt=3 -> exact 66.66…% -> (200 + 1)/3 = 67
                    //  - gc=3, acgt=3 -> exact 100%     -> (300 + 1)/3 = 100 (then clamped to ≤100 below)
                    let gc_percent_bin = ((gc_count as u64 * 100 + (acgt_count as u64 / 2))
                        / acgt_count as u64)
                        .min(100) as usize;

                    counts_by_bin[win_idx].incr(frag_len, gc_percent_bin);
                }
            }

            j += 1;
        }
    }
}
