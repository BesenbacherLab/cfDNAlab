use anyhow::{Result, bail};
use std::sync::Arc;

pub(crate) use crate::shared::length_axis::LengthAxis;

/// Fragment length counter over a shared output axis.
///
/// A `LengthCounts` value stores one count vector over the length axis: either
/// the global result, one genomic window, or one grouped-BED group. Stacking
/// these vectors creates the row dimension of `length_counts.tsv.zst`. The axis
/// is shared so counters can be merged safely without copying the lookup table.
#[derive(Debug, Clone)]
pub(crate) struct LengthCounts {
    pub(crate) counts: Vec<f64>,
    pub(crate) axis: Arc<LengthAxis>,
}

impl LengthCounts {
    /// Create a new zero-initialized `LengthCounts` over a resolved length axis.
    ///
    /// Parameters
    /// ----------
    /// - `axis`:
    ///   Shared output axis used to interpret the count columns.
    pub(crate) fn new(axis: Arc<LengthAxis>) -> Self {
        let counts = vec![0f64; axis.num_bins()];
        Self { counts, axis }
    }

    /// Return the output column for an absolute fragment length.
    ///
    /// Returns `None` when the length is outside the configured axis.
    #[inline]
    pub(crate) fn index_of(&self, length: usize) -> Option<usize> {
        self.axis.bin_index(length)
    }

    /// Increment the counter by a weight for a given fragment length.
    ///
    /// Parameters
    /// ----------
    /// - `length`:
    ///   Absolute fragment length, not a zero-based column index.
    /// - `weight`:
    ///   Count mass to add to the bin containing `length`.
    pub(crate) fn incr_weighted(&mut self, length: usize, weight: f64) -> Result<()> {
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
    #[allow(dead_code)]
    pub(crate) fn get(&self, length: usize) -> Option<f64> {
        self.index_of(length).map(|index| self.counts[index])
    }

    /// Number of length bins.
    #[inline]
    pub(crate) fn n_lengths(&self) -> usize {
        self.counts.len()
    }

    /// Create a zero-initialized `LengthCounts` with the same axis and shape.
    #[inline]
    pub(crate) fn zeroed_like(&self) -> Self {
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
    pub(crate) fn is_compatible_with(&self, other: &Self) -> bool {
        self.counts.len() == other.counts.len()
            && (Arc::ptr_eq(&self.axis, &other.axis) || self.axis.edges() == other.axis.edges())
    }

    /// Merge counts from `other` into `self`.
    ///
    /// The compatibility check guarantees that both counters use the same
    /// column layout before values are added column by column.
    pub(crate) fn merge_from(&mut self, other: &Self) -> Result<()> {
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
    pub(crate) fn collapse<'a, I>(iter: I) -> Result<Self>
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

#[cfg(test)]
mod tests {
    include!("counting_tests.rs");
}
