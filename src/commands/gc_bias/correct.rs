use crate::commands::gc_bias::{
    binning::{BinnedAxis, bins_from_edges},
    counting::{GCPrefixes, get_gc_integer_percentage_for_window},
    package::GCCorrectionPackage,
};
use crate::shared::gc_tag::{SanitizedGCWeight, sanitize_gc_weight};
use crate::shared::interval::Interval;
use crate::shared::reference::twobit_contig_signature;
use anyhow::{Context, Result, anyhow, ensure};
use ndarray::{Array1, Array2, Axis};
use std::str::FromStr;
use tracing::warn;

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
        if self.length_offset_index(fragment_length).is_none() {
            return Ok(None);
        }

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
        Ok(
            match sanitize_gc_weight(self.get_correction_weight(fragment_length, gc_bin)?) {
                SanitizedGCWeight::Usable(weight) => Some(weight),
                SanitizedGCWeight::Unusable { .. } => None,
            },
        )
    }

    /// Get the GC correction weight for a combination of fragment length and GC percentage.
    ///
    /// NOTE: The weight is **multiplicative**, so to correct a fragment's contribution,
    /// **multiply** its existing weight (e.g. `1.0`) with the correction weight.
    #[inline]
    pub fn get_correction_weight(&self, fragment_length: usize, gc_pct: usize) -> Result<f64> {
        let length_bin = self.length_bin(fragment_length)?;
        let gc_bin = self.gc_bin(gc_pct)?;
        Ok(self.correction_matrix[(length_bin, gc_bin)])
    }

    /// Check whether this package covers a fragment length.
    ///
    /// Commands can use this to classify an out-of-package length as an
    /// unusable GC weight before trying to look up a correction.
    ///
    /// Parameters
    /// ----------
    /// - `fragment_length`:
    ///   Aligned fragment length used for GC correction
    ///
    /// Returns
    /// -------
    /// - `bool`:
    ///   `true` when the length is inside the package range
    #[inline]
    pub fn covers_fragment_length(&self, fragment_length: usize) -> bool {
        self.length_offset_index(fragment_length).is_some()
    }

    /// Return the inclusive fragment-length range covered by this package.
    ///
    /// Returns
    /// -------
    /// - `(usize, usize)`:
    ///   Minimum and maximum aligned fragment lengths covered by the package
    #[inline]
    pub fn length_range(&self) -> (usize, usize) {
        (self.length_min, self.length_max)
    }

    #[inline]
    fn length_offset_index(&self, fragment_length: usize) -> Option<usize> {
        if fragment_length < self.length_min || fragment_length > self.length_max {
            return None;
        }
        Some(fragment_length - self.length_min)
    }

    #[inline]
    fn gc_offset_index(&self, gc_pct: usize) -> Option<usize> {
        if gc_pct < self.gc_min || gc_pct > self.gc_max {
            return None;
        }
        Some(gc_pct - self.gc_min)
    }

    #[inline]
    fn length_bin(&self, fragment_length: usize) -> Result<usize> {
        let length_idx = self.length_offset_index(fragment_length).ok_or_else(|| {
            anyhow!(
                "GC correction: unexpected fragment length {}",
                fragment_length
            )
        })?;
        self.lengths_bins
            .index_to_bin
            .get(&length_idx)
            .copied()
            .ok_or_else(|| {
                anyhow!(
                    "GC correction: unexpected fragment length {}",
                    fragment_length
                )
            })
            .with_context(|| format!("length range [{}-{}]", self.length_min, self.length_max))
    }

    #[inline]
    fn gc_bin(&self, gc_pct: usize) -> Result<usize> {
        let gc_idx = self
            .gc_offset_index(gc_pct)
            .ok_or_else(|| anyhow!("GC correction: unexpected GC percentage {}", gc_pct))?;
        self.gc_bins
            .index_to_bin
            .get(&gc_idx)
            .copied()
            .ok_or_else(|| anyhow!("GC correction: unexpected GC percentage {}", gc_pct))
            .with_context(|| format!("GC range [{}-{}]", self.gc_min, self.gc_max))
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
    /// The GC package stores binned length correction curves. Requested-range selection keeps
    /// every package length bin that overlaps the requested inclusive fragment-length range.
    pub fn from_gc_corrector(
        corrector: &GCCorrector,
        weighting_scheme: &MarginalizeLengthsWeightingScheme,
        gc_length_range: GCLengthRange,
        min_fragment_length: u32,
        max_fragment_length: u32,
    ) -> Result<Self> {
        let selected_length_bins = selected_length_bins(
            corrector,
            gc_length_range,
            min_fragment_length,
            max_fragment_length,
        )?;
        let selected_correction_matrix = corrector
            .correction_matrix
            .select(Axis(0), &selected_length_bins);

        // Average corrections to remove the length dimension
        // using the specified weighting scheme
        let correction_vector = match weighting_scheme {
            MarginalizeLengthsWeightingScheme::Equal => selected_correction_matrix
                .mean_axis(Axis(0))
                .ok_or_else(|| anyhow!("No GC correction length bins selected"))?,
            MarginalizeLengthsWeightingScheme::Frequency => {
                let weights = corrector
                    .length_bin_frequencies
                    .select(Axis(0), &selected_length_bins);
                let total_weight: f64 = weights.iter().sum();
                ensure!(total_weight > 0.0, "Length-bin frequencies sum to zero");
                // correction_matrix shape: (length_bins, gc_bins)
                selected_correction_matrix
                    .t()
                    .dot(&weights)
                    .mapv(|v| v / total_weight)
            }
            MarginalizeLengthsWeightingScheme::MaxFrequency => {
                let weights = corrector
                    .length_bin_frequencies
                    .select(Axis(0), &selected_length_bins);
                let (max_index, max_frequency) = weights.iter().copied().enumerate().fold(
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
                selected_correction_matrix
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

fn selected_length_bins(
    corrector: &GCCorrector,
    gc_length_range: GCLengthRange,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<Vec<usize>> {
    match gc_length_range {
        GCLengthRange::Package => Ok((0..corrector.correction_matrix.dim().0).collect()),
        GCLengthRange::Requested => {
            let requested_min = min_fragment_length as usize;
            let requested_max = max_fragment_length as usize;
            ensure!(
                requested_min <= requested_max,
                "minimum fragment length ({}) must be <= maximum fragment length ({})",
                requested_min,
                requested_max
            );
            ensure!(
                requested_min >= corrector.length_min && requested_max <= corrector.length_max,
                "requested GC length range [{}-{}] is outside package range [{}-{}]",
                requested_min,
                requested_max,
                corrector.length_min,
                corrector.length_max
            );

            let first_selected_bin = corrector.length_bin(requested_min)?;
            let last_selected_bin = corrector.length_bin(requested_max)?;

            // Package length bins are contiguous, so every bin between the endpoint bins overlaps
            // the requested inclusive fragment-length range
            let selected_bins: Vec<usize> = (first_selected_bin..=last_selected_bin).collect();
            ensure!(
                !selected_bins.is_empty(),
                "requested GC length range [{}-{}] selected no correction length bins",
                requested_min,
                requested_max
            );
            Ok(selected_bins)
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub enum GCLengthRange {
    #[default]
    Requested,
    Package,
}

impl FromStr for GCLengthRange {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "requested" {
            Ok(GCLengthRange::Requested)
        } else if s == "package" {
            Ok(GCLengthRange::Package)
        } else {
            Err("Use 'requested' or 'package'".into())
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub enum MarginalizeLengthsWeightingScheme {
    #[default]
    Equal,
    Frequency,
    MaxFrequency,
}

impl FromStr for MarginalizeLengthsWeightingScheme {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "equal" {
            Ok(MarginalizeLengthsWeightingScheme::Equal)
        } else if s == "frequency" {
            Ok(MarginalizeLengthsWeightingScheme::Frequency)
        } else if s == "max-frequency" {
            Ok(MarginalizeLengthsWeightingScheme::MaxFrequency)
        } else {
            Err("Use 'equal', 'frequency', or 'max-frequency'".into())
        }
    }
}

pub fn load_gc_corrector<P: AsRef<std::path::Path>>(
    gc_file: Option<&P>,
    ref_2bit: Option<&P>,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<Option<GCCorrector>> {
    if let Some(path) = gc_file {
        let package = GCCorrectionPackage::from_file(path)?;
        warn_on_reference_contig_mismatch(&package, ref_2bit)?;
        validate_gc_package_compatibility(&package, min_fragment_length, max_fragment_length)?;
        Ok(Some(GCCorrector::from_package(&package)?))
    } else {
        Ok(None)
    }
}

pub fn load_length_agnostic_gc_corrector<P: AsRef<std::path::Path>>(
    gc_file: Option<&P>,
    ref_2bit: Option<&P>,
    weighting_scheme: &MarginalizeLengthsWeightingScheme,
    gc_length_range: GCLengthRange,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<Option<LengthAgnosticGCCorrector>> {
    if let Some(path) = gc_file {
        let package = GCCorrectionPackage::from_file(path)?;
        warn_on_reference_contig_mismatch(&package, ref_2bit)?;
        validate_gc_package_compatibility(&package, min_fragment_length, max_fragment_length)?;
        let gc_corrector = GCCorrector::from_package(&package)?;
        let length_agnostic_gc_corrector = LengthAgnosticGCCorrector::from_gc_corrector(
            &gc_corrector,
            weighting_scheme,
            gc_length_range,
            min_fragment_length,
            max_fragment_length,
        )?;
        Ok(Some(length_agnostic_gc_corrector))
    } else {
        Ok(None)
    }
}

fn warn_on_reference_contig_mismatch<P: AsRef<std::path::Path>>(
    package: &GCCorrectionPackage,
    ref_2bit: Option<&P>,
) -> Result<()> {
    let Some(ref_2bit) = ref_2bit else {
        return Ok(());
    };

    let run_signature = twobit_contig_signature(ref_2bit)?;
    if run_signature != package.reference_contig_signature {
        warn!(
            target: "gc-bias",
            "GC correction package was built against a reference contig set that differs from --ref-2bit. GC weights may be invalid if the package and input reference are from different assemblies."
        );
    }
    Ok(())
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

#[cfg(test)]
mod tests {
    include!("correct_tests.rs");
}
