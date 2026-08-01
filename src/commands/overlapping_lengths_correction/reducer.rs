//! Additive statistics that connect parallel tile scanning to the two in-memory model fits.

use std::collections::BTreeMap;

use anyhow::{Result, ensure};

/// Additive sufficient statistics collected by a tile.
///
/// Keeping only these values is what allows both model fits to run after a single tiled sweep of
/// the BAM. The second fit divides each bin's signal sum by the first two bin-wise correction
/// factors. It does not need positional coverage again.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct OverlappingLengthStatistics {
    /// Monotonically increasing bin boundaries, including configured minimum and maximum.
    pub(crate) length_bin_edges: Vec<f64>,
    /// Number of eligible covered genomic bases assigned to each length bin.
    pub(crate) length_bin_base_counts: Vec<u64>,
    /// Sum of observed, optionally corrected coverage signal within each length bin.
    pub(crate) observed_signal_sums: Vec<f64>,
    /// Frequency of each positive raw integer fragment depth over eligible bases.
    pub(crate) raw_depth_frequencies: BTreeMap<u32, u64>,
    /// Total covered, non-blacklisted positions accumulated across all bins.
    pub(crate) eligible_covered_bases: u64,
}

impl OverlappingLengthStatistics {
    /// Allocate an empty reducer for a fixed set of length-bin edges.
    pub(crate) fn new(length_bin_edges: Vec<f64>) -> Result<Self> {
        ensure!(
            length_bin_edges.len() >= 2,
            "average overlapping fragment length bins require at least two edges"
        );
        ensure!(
            length_bin_edges.windows(2).all(|pair| pair[0] < pair[1]),
            "average overlapping fragment length bin edges must increase strictly"
        );
        let num_bins = length_bin_edges.len() - 1;
        Ok(Self {
            length_bin_edges,
            length_bin_base_counts: vec![0; num_bins],
            observed_signal_sums: vec![0.0; num_bins],
            raw_depth_frequencies: BTreeMap::new(),
            eligible_covered_bases: 0,
        })
    }

    /// Add a covered, eligible genomic position to all sufficient statistics.
    ///
    /// `average_overlapping_length` selects the length bin. `observed_signal` contributes to the
    /// curve being corrected. `raw_depth` contributes independently to the sampling-variance
    /// mixture used by both fits.
    pub(crate) fn add_position(
        &mut self,
        average_overlapping_length: f64,
        observed_signal: f64,
        raw_depth: u32,
    ) -> Result<()> {
        ensure!(
            average_overlapping_length.is_finite() && average_overlapping_length > 0.0,
            "average overlapping fragment length must be finite and positive"
        );
        ensure!(
            observed_signal.is_finite() && observed_signal >= 0.0,
            "observed coverage signal must be finite and non-negative"
        );
        ensure!(
            raw_depth > 0,
            "raw depth must be positive for a covered base"
        );

        let bin_index = self.bin_index(average_overlapping_length);
        self.length_bin_base_counts[bin_index] += 1;
        self.observed_signal_sums[bin_index] += observed_signal;
        *self.raw_depth_frequencies.entry(raw_depth).or_default() += 1;
        self.eligible_covered_bases += 1;
        Ok(())
    }

    /// Merge another tile's statistics into this reducer.
    ///
    /// Every stored quantity is additive. Equal bin edges are required so vector positions retain
    /// the same scientific meaning after parallel reduction.
    pub(crate) fn merge(&mut self, other: Self) -> Result<()> {
        ensure!(
            self.length_bin_edges == other.length_bin_edges,
            "cannot merge overlapping-length statistics with different bins"
        );
        for (target, value) in self
            .length_bin_base_counts
            .iter_mut()
            .zip(other.length_bin_base_counts)
        {
            *target += value;
        }
        for (target, value) in self
            .observed_signal_sums
            .iter_mut()
            .zip(other.observed_signal_sums)
        {
            *target += value;
        }
        for (depth, count) in other.raw_depth_frequencies {
            *self.raw_depth_frequencies.entry(depth).or_default() += count;
        }
        self.eligible_covered_bases += other.eligible_covered_bases;
        Ok(())
    }

    /// Calculate the observed mean signal in every configured length bin.
    ///
    /// Empty bins are rejected because the current exact fitting sequence requires a finite value
    /// at every retained bin. This first version deliberately keeps all configured bins instead of
    /// silently dropping extremes.
    pub(crate) fn observed_means(&self) -> Result<Vec<f64>> {
        self.length_bin_base_counts
            .iter()
            .zip(&self.observed_signal_sums)
            .enumerate()
            .map(|(bin_index, (&count, &sum))| {
                ensure!(
                    count > 0,
                    "average overlapping fragment length bin {} has no eligible covered bases",
                    bin_index
                );
                let mean = sum / count as f64;
                ensure!(
                    mean.is_finite(),
                    "average overlapping fragment length bin {} has a non-finite coverage mean",
                    bin_index
                );
                Ok(mean)
            })
            .collect()
    }

    /// Map a positional average length to a bin, clipping values to the extreme bins.
    ///
    /// The configured maximum is included in the final bin. Clipping also mirrors the lookup used
    /// when the fitted package is later applied by fcoverage.
    fn bin_index(&self, value: f64) -> usize {
        let num_bins = self.length_bin_edges.len() - 1;
        let insertion = self.length_bin_edges.partition_point(|edge| *edge <= value);
        insertion.saturating_sub(1).min(num_bins - 1)
    }
}

/// Build length-bin edges anchored at the configured minimum.
///
/// The maximum is always appended as the final edge. Consequently, the final bin is shorter when
/// the inclusive fitting range is not divisible by `width`, and a value equal to `maximum` is
/// assigned to that final bin by `OverlappingLengthStatistics::bin_index`.
pub(crate) fn build_length_bin_edges(minimum: u32, maximum: u32, width: u32) -> Result<Vec<f64>> {
    ensure!(
        minimum <= maximum,
        "minimum fragment length must not exceed maximum"
    );
    ensure!(width > 0, "length bin size must be positive");
    ensure!(
        minimum < maximum,
        "fragment length fitting interval must contain at least two lengths"
    );
    let mut edges = vec![minimum as f64];
    let mut next = minimum.saturating_add(width);
    while next < maximum {
        edges.push(next as f64);
        next = next.saturating_add(width);
    }
    edges.push(maximum as f64);
    Ok(edges)
}

#[cfg(test)]
mod tests {
    include!("reducer_tests.rs");
}
