use anyhow::{Result, ensure};

/// Histogram specification with validated edges and counts.
///
/// Ensures edges are strictly increasing and aligned to the counts length so
/// callers can draw histograms without re-validating inputs. Supports
/// construction from pre-binned counts or raw data with a fixed bin width.
///
/// # Parameters
/// - `edges`:
///     Bin boundaries in ascending order, inclusive of the left edge and
///     exclusive of the right edge for every bin except the final edge which
///     closes the interval.
/// - `counts`:
///     Mass assigned to each interval between successive edges.
///
/// # Returns
/// - `HistogramSpec`:
///     Validated histogram ready for plotting.
#[derive(Clone, Debug)]
pub struct HistogramSpec {
    pub edges: Vec<f64>,
    pub counts: Vec<f64>,
}

impl HistogramSpec {
    /// Build a histogram from counts paired with explicit edges.
    ///
    /// Use this when you already have a count per interval and you want to
    /// supply the exact interval boundaries. Validates that edge and count
    /// lengths agree, edges are strictly increasing, and counts are finite and
    /// non-negative.
    pub fn from_binned(edges: Vec<f64>, counts: Vec<f64>) -> Result<Self> {
        ensure!(
            edges.len() >= 2,
            "Histogram edges must contain at least two entries"
        );
        ensure!(
            edges.len() == counts.len() + 1,
            "Histogram edges must be one longer than counts"
        );
        validate_edges(&edges)?;
        ensure!(
            counts.iter().all(|v| v.is_finite() && *v >= 0.0),
            "Histogram counts must be finite and non-negative"
        );

        Ok(Self { edges, counts })
    }

    /// Build a histogram by binning raw observations into fixed-width bins.
    ///
    /// Uses half-open intervals `[edge_i, edge_{i+1})` with the final bin
    /// including the right edge to avoid dropping values at the boundary.
    /// Rejects values that are non-finite or outside the configured range so
    /// callers notice data quality problems instead of silently discarding
    /// input. Rarely used in this codebase because most callers operate on
    /// pre-summarized counts rather than raw observations.
    pub fn from_data(data: &[f64], start: f64, end: f64, bin_width: f64) -> Result<Self> {
        ensure!(bin_width > 0.0, "Histogram bin width must be positive");
        ensure!(end > start, "Histogram end must be greater than start");

        let mut edges = vec![start];
        let mut current = start;
        while current < end {
            let next = (current + bin_width).min(end);
            ensure!(
                next > current,
                "Bin width produced a non-increasing edge, check bin width and range"
            );
            edges.push(next);
            current = next;
        }

        let mut counts = vec![0.0; edges.len() - 1];
        let mut rejected = 0usize;

        for &value in data.iter() {
            if !value.is_finite() || value < start || value > end {
                rejected += 1;
                continue;
            }
            let mut bin_idx = ((value - start) / bin_width).floor() as usize;
            if bin_idx >= counts.len() {
                bin_idx = counts.len() - 1;
            }
            counts[bin_idx] += 1.0;
        }

        ensure!(
            rejected == 0,
            "Histogram input contained {} invalid values (non-finite or out of range)",
            rejected
        );

        Self::from_binned(edges, counts)
    }

    /// Total mass stored in the histogram.
    pub fn total(&self) -> f64 {
        self.counts.iter().copied().sum()
    }

    /// Build a histogram when you already have counts per equal-width bin.
    ///
    /// Useful for precomputed vectors such as sums across a matrix axis where
    /// bins are uniform and contiguous. When `bin_width` is `None`, a unit
    /// width of 1.0 is used. Prefer this over `from_binned` when you only have
    /// counts and a uniform spacing, not arbitrary edges.
    pub fn from_counts(counts: Vec<f64>, bin_width: Option<f64>) -> Result<Self> {
        ensure!(
            !counts.is_empty(),
            "Histogram counts must contain at least one entry"
        );
        let width = bin_width.unwrap_or(1.0);
        ensure!(width > 0.0, "Histogram bin width must be positive");
        let edges: Vec<f64> = (0..=counts.len())
            .map(|i| i as f64 * width)
            .collect();
        Self::from_binned(edges, counts)
    }

    /// Maximum bin mass.
    pub fn max(&self) -> f64 {
        self.counts.iter().copied().fold(0.0, |acc, v| acc.max(v))
    }
}

fn validate_edges(edges: &[f64]) -> Result<()> {
    for window in edges.windows(2) {
        let left = window[0];
        let right = window[1];
        ensure!(
            right > left,
            "Histogram edges must be strictly increasing: found {} then {}",
            left,
            right
        );
    }
    Ok(())
}
