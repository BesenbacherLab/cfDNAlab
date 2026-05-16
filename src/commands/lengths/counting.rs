use anyhow::{Result, bail, ensure};
use ndarray::Array2;
use std::sync::Arc;

pub use crate::shared::length_axis::LengthAxis;

/// Fragment-length counter over a shared output axis.
///
/// A `LengthCounts` value stores one count vector over the length axis: either
/// the global result, one genomic window, or one grouped-BED group. Stacking
/// these vectors creates the row dimension of `length_counts.tsv.zst`. The axis
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
