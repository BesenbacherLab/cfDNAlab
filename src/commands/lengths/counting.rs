use crate::commands::cli_common::{MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION};
use anyhow::{Result, bail, ensure};
use ndarray::Array2;
use std::sync::Arc;

/// Resolved fragment length output axis for `cfdna lengths`.
///
/// `lengths` stores output columns as half-open length bins. The axis owns
/// those edges and a lookup table from absolute fragment length to output
/// column. Exact per-bp bins and wider bins therefore use the same column
/// resolution path.
///
/// The lookup table is intentionally shared through `Arc<LengthAxis>` by
/// `LengthCounts`. A run may allocate many thousands of counters, so the lookup
/// table is stored once while each counter owns only its per-bin counts and a
/// pointer to the shared axis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LengthAxis {
    edges: Vec<u32>,
    length_to_bin: Vec<usize>,
    single_bp_bins: bool,
}

impl LengthAxis {
    /// Create a validated length axis from half-open bin edges.
    ///
    /// The first edge is the minimum included fragment length. The final edge
    /// is exclusive, so the maximum included fragment length is `last_edge - 1`.
    /// Every intermediate edge must be strictly increasing.
    ///
    /// Parameters
    /// ----------
    /// - `edges`:
    ///   Half-open bin edges, e.g. `[30, 40, 50]` for `[30,40)` and `[40,50)`.
    ///
    /// Returns
    /// -------
    /// - `LengthAxis`:
    ///   A resolved axis with O(1) length-to-column lookup.
    pub fn new(edges: Vec<u32>) -> Result<Self> {
        ensure!(
            edges.len() >= 2,
            "length bin edges must contain at least two values"
        );
        ensure!(
            edges[0] >= MIN_ACGT_BASES_FOR_GC_FRACTION,
            "length bin edges must be >= {}",
            MIN_ACGT_BASES_FOR_GC_FRACTION
        );
        ensure!(
            *edges.last().expect("length edges checked non-empty")
                <= MAX_SUPPORTED_FRAGMENT_LENGTH + 1,
            "length bin edges must be <= {}",
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        );
        ensure!(
            edges.windows(2).all(|window| window[0] < window[1]),
            "length bin edges must be strictly increasing"
        );

        let single_bp_bins = edges.windows(2).all(|window| window[1] == window[0] + 1);
        let num_bins = edges.len() - 1;
        let max_edge = *edges.last().expect("length edges checked non-empty") as usize;
        let mut length_to_bin = vec![usize::MAX; max_edge.max(1)];

        for bin_index in 0..num_bins {
            let bin_start = edges[bin_index] as usize;
            let bin_end = edges[bin_index + 1] as usize;
            for length in bin_start..bin_end {
                length_to_bin[length] = bin_index;
            }
        }

        Ok(Self {
            edges,
            length_to_bin,
            single_bp_bins,
        })
    }

    /// Return the half-open bin edges.
    ///
    /// Consecutive edge pairs define output columns. The interval for column
    /// `i` is `[edges[i], edges[i + 1])`.
    #[inline]
    pub fn edges(&self) -> &[u32] {
        &self.edges
    }

    /// Number of output length bins.
    #[inline]
    pub fn num_bins(&self) -> usize {
        self.edges.len() - 1
    }

    /// Minimum included fragment length.
    #[inline]
    pub fn min_fragment_length(&self) -> u32 {
        self.edges[0]
    }

    /// Maximum included fragment length.
    #[inline]
    pub fn max_fragment_length(&self) -> u32 {
        *self.edges.last().expect("length edges checked non-empty") - 1
    }

    /// Check whether an absolute fragment length belongs to this output axis.
    ///
    /// This uses the same lookup table as `bin_index()` and is intended for
    /// cheap inclusion checks before adding count mass.
    #[inline]
    pub fn contains(&self, length: u32) -> bool {
        self.bin_index(length as usize).is_some()
    }

    /// Return the output column for an absolute fragment length.
    ///
    /// Returns `None` for lengths below the first edge, at or beyond the final
    /// exclusive edge, or in gaps. Gaps are not produced by the current resolver,
    /// but the sentinel keeps the lookup robust if construction changes later.
    #[inline]
    pub fn bin_index(&self, length: usize) -> Option<usize> {
        let bin_index = *self.length_to_bin.get(length)?;
        (bin_index != usize::MAX).then_some(bin_index)
    }

    /// Whether every length bin has width 1 bp.
    #[inline]
    pub fn is_single_bp_bins(&self) -> bool {
        self.single_bp_bins
    }
}

/// Fragment-length counter over a shared output axis.
///
/// A `LengthCounts` value stores one count vector over the length axis: either
/// the global result, one genomic window, or one grouped-BED group. Stacking
/// these vectors creates the first dimension of `length_counts.npy`. The axis
/// is shared so counters can be merged safely without copying the lookup table.
#[derive(Debug, Clone)]
pub struct LengthCounts {
    pub counts: Vec<f64>,
    pub axis: Arc<LengthAxis>,
}

impl LengthCounts {
    /// Create a new zero-initialized `LengthCounts` over a resolved length axis.
    ///
    /// Parameters
    /// ----------
    /// - `axis`:
    ///   Shared output axis used to interpret the count columns.
    pub fn new(axis: Arc<LengthAxis>) -> Self {
        let counts = vec![0f64; axis.num_bins()];
        Self { counts, axis }
    }

    /// Return the output column for an absolute fragment length.
    ///
    /// Returns `None` when the length is outside the configured axis.
    #[inline]
    pub fn index_of(&self, length: usize) -> Option<usize> {
        self.axis.bin_index(length)
    }

    /// Increment the counter by `1.0` for a given fragment length.
    ///
    /// Returns an error when `length` does not map to any configured bin. The
    /// caller is responsible for using the same axis for filtering and counting.
    pub fn incr(&mut self, length: usize) -> Result<()> {
        self.incr_weighted(length, 1.0)
    }

    /// Increment the counter by a weight for a given fragment length.
    ///
    /// Parameters
    /// ----------
    /// - `length`:
    ///   Absolute fragment length, not a zero-based column index.
    /// - `weight`:
    ///   Count mass to add to the bin containing `length`.
    pub fn incr_weighted(&mut self, length: usize, weight: f64) -> Result<()> {
        let Some(index) = self.index_of(length) else {
            bail!("fragment length {length} did not map to any configured length bin");
        };
        self.counts[index] += weight;
        Ok(())
    }

    /// Get the count for the bin that contains an absolute fragment length.
    ///
    /// With wider bins, multiple lengths return the same value. For example,
    /// with edges `[30,40,50]`, both `get(35)` and `get(39)` read column 0.
    pub fn get(&self, length: usize) -> Option<f64> {
        self.index_of(length).map(|index| self.counts[index])
    }

    /// Number of length bins.
    #[inline]
    pub fn n_lengths(&self) -> usize {
        self.counts.len()
    }

    /// Create a zero-initialized `LengthCounts` with the same axis and shape.
    #[inline]
    pub fn zeroed_like(&self) -> Self {
        Self {
            counts: vec![0f64; self.n_lengths()],
            axis: Arc::clone(&self.axis),
        }
    }

    /// Check whether two counters can be merged without changing column meaning.
    ///
    /// Compatibility requires both the count width and the full axis definition
    /// to match. Same-width axes such as `[10,20]` and `[20,30]` must not merge.
    /// Counters produced inside one run usually share the same `Arc`, so pointer
    /// equality handles the common case without comparing edge vectors.
    #[inline]
    pub fn is_compatible_with(&self, other: &Self) -> bool {
        self.counts.len() == other.counts.len()
            && (Arc::ptr_eq(&self.axis, &other.axis) || self.axis.edges() == other.axis.edges())
    }

    /// Merge counts from `other` into `self`.
    ///
    /// The compatibility check guarantees that both counters use the same
    /// column layout before values are added column by column.
    pub fn merge_from(&mut self, other: &Self) -> Result<()> {
        if !self.is_compatible_with(other) {
            bail!(
                "incompatible LengthCounts: self={} vs other={}",
                self,
                other
            );
        }
        for (index, count_other) in other.counts.iter().enumerate() {
            self.counts[index] += *count_other;
        }
        Ok(())
    }

    /// Collapse compatible counters into one summed counter.
    ///
    /// Returns an error for empty input or incompatible axes.
    pub fn collapse<'a, I>(iter: I) -> Result<Self>
    where
        I: IntoIterator<Item = &'a LengthCounts>,
    {
        let mut iterator = iter.into_iter();
        let first = iterator
            .next()
            .ok_or_else(|| anyhow::anyhow!("collapse requires at least one LengthCounts"))?;
        let mut accumulator = first.clone();
        for counts in iterator {
            accumulator.merge_from(counts)?;
        }
        Ok(accumulator)
    }
}

impl Default for LengthCounts {
    fn default() -> Self {
        let default_edges: Vec<u32> = (30..=1001).collect();
        let axis = Arc::new(LengthAxis::new(default_edges).expect("default length axis is valid"));
        Self::new(axis)
    }
}

impl std::fmt::Display for LengthCounts {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "LengthCounts(lengths:[{}..={}], dims:({}) )",
            self.axis.min_fragment_length(),
            self.axis.max_fragment_length(),
            self.n_lengths(),
        )
    }
}

/// Stack count vectors into a 2D count matrix.
///
/// Each input counter becomes one vector along the first matrix dimension.
/// Every counter must use the same axis, which prevents producing a matrix
/// with identical widths but incompatible length-bin definitions.
pub fn stack_length_counts(all_counts: &[LengthCounts]) -> Result<Array2<f64>> {
    let first = all_counts
        .first()
        .ok_or_else(|| anyhow::anyhow!("stack_length_counts requires at least one counter"))?;
    let num_rows = all_counts.len();
    let num_columns = first.counts.len();

    let mut array = Array2::<f64>::zeros((num_rows, num_columns));
    for (row_index, length_counts) in all_counts.iter().enumerate() {
        ensure!(
            length_counts.is_compatible_with(first),
            "length count entry {} has incompatible length axis",
            row_index
        );
        for (column_index, &count) in length_counts.counts.iter().enumerate() {
            array[(row_index, column_index)] = count;
        }
    }

    Ok(array)
}

#[cfg(test)]
mod tests {
    include!("counting_tests.rs");
}
