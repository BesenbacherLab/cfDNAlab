use crate::shared::constants::{MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION};
use anyhow::{Result, ensure};
use serde::Serialize;

#[derive(Serialize)]
pub(crate) struct LengthAxisSettings<'a> {
    column_intervals: &'static str,
    min_fragment_length: u32,
    max_fragment_length: u32,
    n_bins: usize,
    single_bp_bins: bool,
    bin_definition: LengthBinDefinition<'a>,
}

#[derive(Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum LengthBinDefinition<'a> {
    SteppedRange { start: u32, end: u32, step: u32 },
    ExplicitEdges { edges: &'a [u32] },
}

/// Resolved fragment length output axis.
///
/// Commands store fragment length output columns as half-open length bins. The axis owns those
/// edges and a lookup table from absolute fragment length to output column. Exact per-bp bins and
/// wider bins therefore use the same column resolution path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LengthAxis {
    edges: Vec<u32>,
    length_to_bin: Vec<usize>,
    single_bp_bins: bool,
}

impl LengthAxis {
    /// Create a validated length axis from half-open bin edges.
    ///
    /// The first edge is the minimum included fragment length. The final edge is exclusive, so the
    /// maximum included fragment length is `last_edge - 1`. Every intermediate edge must be
    /// strictly increasing.
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
    pub(crate) fn new(edges: Vec<u32>) -> Result<Self> {
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
            edges
                .windows(2)
                .all(|edge_window| edge_window[0] < edge_window[1]),
            "length bin edges must be strictly increasing"
        );

        let single_bp_bins = edges
            .windows(2)
            .all(|edge_window| edge_window[1] == edge_window[0] + 1);
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
    /// Consecutive edge pairs define output columns. The interval for column `i` is
    /// `[edges[i], edges[i + 1])`.
    #[inline]
    pub(crate) fn edges(&self) -> &[u32] {
        &self.edges
    }

    /// Number of output length bins.
    #[inline]
    pub(crate) fn num_bins(&self) -> usize {
        self.edges.len() - 1
    }

    /// Minimum included fragment length.
    #[inline]
    pub(crate) fn min_fragment_length(&self) -> u32 {
        self.edges[0]
    }

    /// Maximum included fragment length.
    #[inline]
    pub(crate) fn max_fragment_length(&self) -> u32 {
        *self.edges.last().expect("length edges checked non-empty") - 1
    }

    /// Check whether an absolute fragment length belongs to this output axis.
    ///
    /// This uses the same lookup table as `bin_index()` and is intended for cheap inclusion checks
    /// before adding count mass.
    #[inline]
    pub(crate) fn contains(&self, length: u32) -> bool {
        self.bin_index(length as usize).is_some()
    }

    /// Return the output column for an absolute fragment length.
    ///
    /// Returns `None` for lengths below the first edge, at or beyond the final exclusive edge, or
    /// in gaps. Gaps are not produced by the current resolver, but the sentinel keeps the lookup
    /// robust if construction changes later.
    #[inline]
    pub(crate) fn bin_index(&self, length: usize) -> Option<usize> {
        let bin_index = *self.length_to_bin.get(length)?;
        (bin_index != usize::MAX).then_some(bin_index)
    }

    /// Whether every length bin has width 1 bp.
    #[inline]
    pub(crate) fn is_single_bp_bins(&self) -> bool {
        self.single_bp_bins
    }

    /// Return a JSON-serializable description of the length-bin output axis.
    ///
    /// Commands that write dense or table-like length-binned outputs use this shared length-axis
    /// representation so downstream readers see the same half-open bin contract everywhere.
    pub(crate) fn settings(&self) -> LengthAxisSettings<'_> {
        LengthAxisSettings {
            column_intervals: "half_open",
            min_fragment_length: self.min_fragment_length(),
            max_fragment_length: self.max_fragment_length(),
            n_bins: self.num_bins(),
            single_bp_bins: self.is_single_bp_bins(),
            bin_definition: length_bin_definition(self.edges()),
        }
    }
}

fn length_bin_definition(edges: &[u32]) -> LengthBinDefinition<'_> {
    if let Some(step) = stepped_range_step(edges) {
        LengthBinDefinition::SteppedRange {
            start: edges[0],
            end: *edges.last().expect("validated length axis has edges"),
            step,
        }
    } else {
        LengthBinDefinition::ExplicitEdges { edges }
    }
}

fn stepped_range_step(edges: &[u32]) -> Option<u32> {
    let step = edges.get(1)?.checked_sub(*edges.first()?)?;
    if step == 0 {
        return None;
    }

    let last_bin_index = edges.len().checked_sub(2)?;
    for (bin_index, edge_pair) in edges.windows(2).enumerate() {
        let width = edge_pair[1].checked_sub(edge_pair[0])?;
        if bin_index == last_bin_index {
            if width == 0 || width > step {
                return None;
            }
        } else if width != step {
            return None;
        }
    }

    Some(step)
}

#[cfg(test)]
mod tests {
    include!("length_axis_tests.rs");
}
