use crate::commands::gc_bias::{
    binning::{BinnedAxis, bins_from_edges},
    counting::{GCPrefixes, get_gc_integer_percentage_for_window},
    package::GCCorrectionPackage,
};
use crate::shared::interval::Interval;
use anyhow::{Context, Result, anyhow, ensure};
use ndarray::{Array1, Array2, Axis};
use std::str::FromStr;

#[derive(Debug, Clone)]
pub struct GCCorrector {
    correction_matrix: Array2<f64>,
    length_bin_frequencies: Array1<f64>,
    lengths_bins: BinnedAxis,
    gc_bins: BinnedAxis,
    length_min: usize,
    length_max: usize,
    gc_min: usize,
    gc_max: usize,
    end_offset: u64,
}

impl GCCorrector {
    /// Create a `GCCorrector` instance from a loaded `GCCorrectionPackage`.
    pub fn from_package(package: &GCCorrectionPackage) -> Result<Self> {
        let length_bins = bins_from_edges(&package.length_edges)?;
        let gc_bins = bins_from_edges(&package.gc_edges)?;
        let length_min = *package
            .length_edges
            .first()
            .ok_or_else(|| anyhow!("GC correction package contained no length edges"))?;
        let length_max = *package
            .length_edges
            .last()
            .ok_or_else(|| anyhow!("GC correction package contained no length edges"))?;
        let gc_min = *package
            .gc_edges
            .first()
            .ok_or_else(|| anyhow!("GC correction package contained no GC edges"))?;
        let gc_max = *package
            .gc_edges
            .last()
            .ok_or_else(|| anyhow!("GC correction package contained no GC edges"))?;

        Ok(GCCorrector {
            correction_matrix: package.correction_matrix.clone(),
            length_bin_frequencies: package.length_bin_frequencies.clone(),
            lengths_bins: length_bins,
            gc_bins,
            length_min: length_min as usize,
            length_max: length_max as usize,
            gc_min: gc_min as usize,
            gc_max: gc_max as usize,
            end_offset: package.end_offset,
        })
    }

    /// Get multiplicative GC correction weight from fragment coordinates and tile-/chromosome-wise prefix arrays
    ///
    /// **NOTE**: Coordinates must be relative to the prefix arrays.
    #[inline]
    pub fn correct_fragment(
        &self,
        fragment_interval: Interval<u64>,
        gc_prefixes: &GCPrefixes,
    ) -> Result<Option<f64>> {
        let fragment_length = fragment_interval.len() as usize;
        let gc_window = fragment_interval
            .contract(self.end_offset)
            .ok_or_else(|| {
                anyhow!(
                    "GC correction: After applying end-offsets the fragment has no bases left to count GCs at.\
                    Does the minimum fragment length match the one in the reference bias?"
                )
            })?
            .try_to_usize()?;
        let gc_bin = if let Some(gc_pct) =
            get_gc_integer_percentage_for_window(gc_prefixes, gc_window, 0.0, 10)?
        {
            gc_pct
        } else {
            return Ok(None);
        };
        Ok(Some(self.get_correction_weight(fragment_length, gc_bin)?))
    }

    /// Get the GC correction weight for a combination of fragment length and GC percentage.
    ///
    /// NOTE: The weight is **multiplicative**, so to correct a fragment's contribution,
    /// **multiply** its existing weight (e.g. `1.0`) with the correction weight.
    #[inline]
    pub fn get_correction_weight(&self, fragment_length: usize, gc_pct: usize) -> Result<f64> {
        // The length index and GC index have the minimum values assigned at 0
        // So we offset the values by their minimum to get the index, which
        // in turn will give us the bin (for which we have correction factors)
        let length_idx = fragment_length - self.length_min;
        let gc_idx = gc_pct - self.gc_min;
        let length_bin = self
            .lengths_bins
            .index_to_bin
            .get(&length_idx)
            .copied()
            .ok_or_else(|| {
                anyhow!(
                    "GC correction: unexpected fragment length {}",
                    fragment_length
                )
            })
            .with_context(|| format!("length range [{}-{}]", self.length_min, self.length_max))?;

        let gc_bin = self
            .gc_bins
            .index_to_bin
            .get(&gc_idx)
            .copied()
            .ok_or_else(|| anyhow!("GC correction: unexpected GC percentage {}", gc_pct))
            .with_context(|| format!("GC range [{}-{}]", self.gc_min, self.gc_max))?;

        Ok(self.correction_matrix[(length_bin, gc_bin)])
    }
}

#[derive(Debug, Clone)]
pub struct LengthAgnosticGCCorrector {
    correction_vector: Array1<f64>,
    gc_bins: BinnedAxis,
    gc_min: usize,
    gc_max: usize,
    end_offset: u64,
}

impl LengthAgnosticGCCorrector {
    /// Create length-agnostic GC corrector from the standard GC corrector.
    ///
    /// Averages out the length dimension, weighted based on the `weighting_scheme`.
    pub fn from_gc_corrector(
        corrector: &GCCorrector,
        weighting_scheme: &MarginalizeLengthsWeightingScheme,
    ) -> Result<Self> {
        // Average corrections to remove the length dimension
        // using the specified weighting scheme
        let correction_vector = match weighting_scheme {
            MarginalizeLengthsWeightingScheme::Equal => corrector
                .correction_matrix
                .mean_axis(Axis(0))
                .ok_or_else(|| {
                    anyhow!(
                        "Failed to average out the length dimension of the GC correction matrix."
                    )
                })?,
            MarginalizeLengthsWeightingScheme::Coverage => {
                let weights = &corrector.length_bin_frequencies;
                let total_weight: f64 = weights.iter().sum();
                ensure!(total_weight > 0.0, "Length-bin frequencies sum to zero");
                // correction_matrix shape: (length_bins, gc_bins)
                corrector
                    .correction_matrix
                    .t()
                    .dot(weights)
                    .mapv(|v| v / total_weight)
            }
            MarginalizeLengthsWeightingScheme::MaxCoverage => {
                let length_bin_frequencies = &corrector.length_bin_frequencies;
                let (max_index, max_frequency) =
                    length_bin_frequencies.iter().copied().enumerate().fold(
                        (None, 0.0_f64),
                        |(best_index, best_frequency), (index, frequency)| {
                            if frequency > best_frequency {
                                (Some(index), frequency)
                            } else {
                                (best_index, best_frequency)
                            }
                        },
                    );
                let most_frequent_index =
                    max_index.ok_or_else(|| anyhow!("Length-bin frequencies array is empty"))?;
                ensure!(max_frequency > 0.0, "Length-bin frequencies sum to zero");
                corrector
                    .correction_matrix
                    .row(most_frequent_index)
                    .to_owned()
            }
        };

        Ok(Self {
            correction_vector,
            gc_bins: corrector.gc_bins.clone(),
            gc_min: corrector.gc_min,
            gc_max: corrector.gc_max,
            end_offset: corrector.end_offset,
        })
    }

    /// Get multiplicative GC correction weight from fragment coordinates and tile-/chromosome-wise prefix arrays
    ///
    /// **NOTE**: Coordinates must be relative to the prefix arrays.
    #[inline]
    pub fn correct_fragment(
        &self,
        fragment_interval: Interval<u64>,
        gc_prefixes: &GCPrefixes,
    ) -> Result<Option<f64>> {
        let gc_window = fragment_interval
            .contract(self.end_offset)
            .ok_or_else(|| {
                anyhow!(
                    "GC correction: After applying end-offsets the fragment has no bases left to count GCs at.\
                    Does the minimum fragment length match the one in the reference bias?"
                )
            })?
            .try_to_usize()?;
        let gc_bin = if let Some(gc_pct) =
            get_gc_integer_percentage_for_window(gc_prefixes, gc_window, 0.0, 10)?
        {
            gc_pct
        } else {
            return Ok(None);
        };
        Ok(Some(self.get_correction_weight(gc_bin)?))
    }

    /// Get the GC correction weight for a GC percentage.
    ///
    /// NOTE: The weight is **multiplicative**, so to correct a fragment's contribution,
    /// **multiply** its existing weight (e.g. `1.0`) with the correction weight.
    #[inline]
    pub fn get_correction_weight(&self, gc_pct: usize) -> Result<f64> {
        // The GC index has the minimum value assigned at 0
        // So we offset the values by their minimum to get the index, which
        // in turn will give us the bin (for which we have correction factors)
        let gc_idx = gc_pct - self.gc_min;
        let gc_bin = self
            .gc_bins
            .index_to_bin
            .get(&gc_idx)
            .copied()
            .ok_or_else(|| anyhow!("GC correction: unexpected GC percentage {}", gc_pct))
            .with_context(|| format!("GC range [{}-{}]", self.gc_min, self.gc_max))?;

        Ok(self.correction_vector[gc_bin])
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub enum MarginalizeLengthsWeightingScheme {
    #[default]
    Equal,
    Coverage,
    MaxCoverage,
}

impl FromStr for MarginalizeLengthsWeightingScheme {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "equal" {
            Ok(MarginalizeLengthsWeightingScheme::Equal)
        } else if s == "coverage" {
            Ok(MarginalizeLengthsWeightingScheme::Coverage)
        } else if s == "max-coverage" {
            Ok(MarginalizeLengthsWeightingScheme::MaxCoverage)
        } else {
            Err("Use 'equal', 'coverage', or 'max-coverage'".into())
        }
    }
}

pub fn load_gc_corrector<P: AsRef<std::path::Path>>(
    gc_file: Option<&P>,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<Option<GCCorrector>> {
    if let Some(path) = gc_file {
        let package = GCCorrectionPackage::from_file(path)?;
        validate_gc_package_compatibility(&package, min_fragment_length, max_fragment_length)?;
        Ok(Some(GCCorrector::from_package(&package)?))
    } else {
        Ok(None)
    }
}

pub fn load_length_agnostic_gc_corrector<P: AsRef<std::path::Path>>(
    gc_file: Option<&P>,
    weighting_scheme: &MarginalizeLengthsWeightingScheme,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<Option<LengthAgnosticGCCorrector>> {
    if let Some(path) = gc_file {
        let package = GCCorrectionPackage::from_file(path)?;
        validate_gc_package_compatibility(&package, min_fragment_length, max_fragment_length)?;
        let gc_corrector = GCCorrector::from_package(&package)?;
        let length_agnostic_gc_corrector =
            LengthAgnosticGCCorrector::from_gc_corrector(&gc_corrector, weighting_scheme)?;
        Ok(Some(length_agnostic_gc_corrector))
    } else {
        Ok(None)
    }
}

fn validate_gc_package_compatibility(
    package: &GCCorrectionPackage,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<()> {
    let package_min_length = *package
        .length_edges
        .first()
        .context("GC correction package contained no length edges")?;
    let package_max_length = *package
        .length_edges
        .last()
        .context("GC correction package contained no length edges")?;

    let min_length = min_fragment_length;
    let max_length = max_fragment_length;

    let end_offset_twice = 2 * package.end_offset as u32;

    ensure!(
        min_length > end_offset_twice,
        "GC correction: minimum fragment length ({min_length}) must exceed twice the end-offset ({}) used when building the correction. \
        Increase --min-fragment-length or rebuild the GC correction package with a smaller --end-offset.",
        end_offset_twice
    );

    ensure!(
        min_length >= package_min_length && max_length <= package_max_length,
        "GC correction: fragment length range [{}-{}] is outside the range covered by the correction package [{}-{}]. \
        Adjust --min-fragment-length/--max-fragment-length or rebuild the GC correction package with matching limits.",
        min_fragment_length,
        max_fragment_length,
        package_min_length,
        package_max_length
    );

    Ok(())
}
