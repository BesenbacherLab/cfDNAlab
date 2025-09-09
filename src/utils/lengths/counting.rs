use ndarray::Array2;

/// Count array for fragment coverage across fragment lengths.
#[derive(Debug, Clone)]
pub struct LengthCounts {
    pub counts: Vec<u64>,
    pub length_min: usize,
    pub length_max: usize,
}

impl LengthCounts {
    /// Create a new `LengthCounts` with specified length range.
    ///
    /// Parameters
    /// ----------
    /// length_min: usize
    ///     Minimum fragment length (inclusive).
    /// length_max: usize
    ///     Maximum fragment length (inclusive).
    ///
    /// Returns
    /// -------
    /// counts: LengthCounts
    ///     A `LengthCounts` object with all counts initialized to zero.
    pub fn new(length_min: usize, length_max: usize) -> Self {
        let num_lengths = length_max - length_min + 1;
        let counts = vec![0u64; num_lengths];
        Self {
            counts,
            length_min,
            length_max,
        }
    }

    /// Check whether `length` is within configured range.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length to test.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     True if length is in range.
    #[inline]
    fn in_bounds(&self, length: usize) -> bool {
        (self.length_min..=self.length_max).contains(&length)
    }

    /// Compute indices for `length` if in bounds.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute, not zero-based).
    ///
    /// Returns
    /// -------
    /// idx: Option<usize>
    ///     Zero-based index if in range, otherwise `None`.
    #[inline]
    pub fn index_of(&self, length: usize) -> Option<usize> {
        if self.in_bounds(length) {
            Some(length - self.length_min)
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
    pub fn incr(&mut self, length: usize) {
        if let Some(i) = self.index_of(length) {
            self.counts[i] = self.counts[i].saturating_add(1);
        }
    }

    // Get the count at a given fragment length bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    ///
    /// Returns
    /// -------
    /// count: Option<u64>
    ///     The count if index is in range, otherwise `None`.
    pub fn get(&self, length: usize) -> Option<u64> {
        self.index_of(length).map(|i| self.counts[i])
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

    /// Create a zero-initialized `LengthCounts` with the same range and shape as `self`.
    ///
    /// Returns
    /// -------
    /// counts: LengthCounts
    ///     New object with all counts set to zero and identical configuration.
    #[inline]
    pub fn zeroed_like(&self) -> Self {
        let n_len = self.n_lengths();
        Self {
            counts: vec![0u64; n_len],
            length_min: self.length_min,
            length_max: self.length_max,
        }
    }

    /// Check if two `LengthCounts` are compatible for merging (same range).
    ///
    /// Parameters
    /// ----------
    /// other: &LengthCounts
    ///     The other counts array to compare with.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     `true` if the two objects have identical range.
    #[inline]
    pub fn is_compatible_with(&self, other: &Self) -> bool {
        self.length_min == other.length_min && self.length_max == other.length_max
    }

    /// Merge (sum) counts from `other` into `self` using saturating addition.
    ///
    /// Parameters
    /// ----------
    /// other: &LengthCounts
    ///     Counts to add into `self`. Must be compatible (same range/shape).
    ///
    /// Returns
    /// -------
    /// result: Result<(), anyhow::Error>
    ///     An error is returned if the two objects are incompatible.
    pub fn merge_from(&mut self, other: &Self) -> anyhow::Result<()> {
        if !self.is_compatible_with(other) {
            return Err(anyhow::anyhow!(
                "incompatible LengthCounts: self={} vs other={}",
                self,
                other
            ));
        }
        for (i, count_other) in other.counts.iter().enumerate() {
            self.counts[i] = self.counts[i].saturating_add(*count_other);
        }
        Ok(())
    }

    /// Collapse (sum) an iterator of `LengthCounts` into a single object.
    ///
    /// All inputs must be compatible (same ranges/shape). Uses saturating addition.
    ///
    /// Parameters
    /// ----------
    /// iter: IntoIterator<Item = &LengthCounts>
    ///     Collection of references to `LengthCounts` to be summed.
    ///
    /// Returns
    /// -------
    /// total: LengthCounts
    ///     The element-wise sum across all inputs.
    ///
    /// Examples
    /// --------
    /// ```ignore
    /// // Sum across chromosomes:
    /// // let by_chr: HashMap<String, LengthCounts> = ...;
    /// let total = LengthCounts::collapse(by_chr.values())?;
    /// ```
    pub fn collapse<'a, I>(iter: I) -> anyhow::Result<Self>
    where
        I: IntoIterator<Item = &'a LengthCounts>,
    {
        let mut it = iter.into_iter();
        let first = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("collapse requires at least one LengthCounts"))?;
        let mut acc = first.clone(); // Start from the first
        for g in it {
            acc.merge_from(g)?; // Then merge the rest
        }
        Ok(acc)
    }
}

impl Default for LengthCounts {
    /// Create an empty default `LengthCounts` (0–100 GC, 20–600 length).
    fn default() -> Self {
        Self::new(20, 600)
    }
}

impl std::fmt::Display for LengthCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "LengthCounts(len:[{}..={}], dims:({}) )",
            self.length_min,
            self.length_max,
            self.n_lengths(),
        )
    }
}

/// Stack counts from vector of `LengthCounts` to a single 2d array.
pub fn stack_length_counts(all_counts: &Vec<LengthCounts>) -> Array2<u64> {
    // Assume all LengthCounts.counts are the same length
    let n = all_counts.len();
    let m = all_counts[0].counts.len();

    // Allocate a 2D array
    let mut arr = Array2::<u64>::zeros((n, m));

    // Fill the array
    for (i, lc) in all_counts.iter().enumerate() {
        for (j, &c) in lc.counts.iter().enumerate() {
            arr[(i, j)] = c;
        }
    }

    arr
}
