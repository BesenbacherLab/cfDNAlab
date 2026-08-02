//! Versioned Zarr representation shared by model fitting and fcoverage application.

use std::{path::Path, sync::Arc};

use anyhow::{Context, Result, bail, ensure};
use serde_json::{Value, json};
use zarrs::array::data_type;
use zarrs::filesystem::FilesystemStore;

use crate::shared::zarr::{
    create_zarr_store, ensure_zarr_schema, read_zarr_array1, read_zarr_root_attributes,
    write_single_chunk_zarr_array, write_zarr_root_metadata,
};

use super::{
    config::OverlappingLengthsCorrectionConfig,
    model::{
        BASE_SIGMA, OverlappingLengthModel, STUDENT_T_DEGREES_OF_FREEDOM,
        TARGET_MEAN_FRAGMENT_LENGTH,
    },
    reducer::OverlappingLengthStatistics,
};

/// Current on-disk schema version for overlapping fragment length model packages.
pub const OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION: u32 = 2;
/// Schema identifier stored in the Zarr root attributes.
const OVERLAPPING_LENGTHS_CORRECTION_SCHEMA: &str = "overlap_length_model";

/// Three free parameters of a fitted depth-mixture distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureFitParameters {
    /// Multiplier applied to the empirical base sigma at every raw depth.
    pub scale_multiplier: f64,
    /// Coefficient used by LIONHEART's error-function skew term.
    pub skewness: f64,
    /// Location parameter of every depth-specific Student-t component.
    pub mean_fragment_length: f64,
}

/// Self-contained average overlapping fragment length normalization model.
///
/// Source paths and file fingerprints are deliberately not persisted. The metadata records only
/// scientific settings that help a later fcoverage run identify potentially inconsistent modes.
#[derive(Clone, Debug, PartialEq)]
pub struct OverlappingLengthsCorrectionPackage {
    /// On-disk schema version.
    pub version: u32,
    /// Average overlapping fragment length bin boundaries.
    pub length_bin_edges: Vec<f64>,
    /// Midpoints derived from adjacent `length_bin_edges`.
    pub length_bin_midpoints: Vec<f64>,
    /// Number of eligible covered bases assigned to each length bin.
    pub length_bin_base_counts: Vec<u64>,
    /// Sum of the observed coverage signal assigned to each length bin.
    pub observed_signal_sums: Vec<f64>,
    /// Sorted positive raw depths represented in the depth histogram.
    pub raw_depths: Vec<u32>,
    /// Genomic base frequencies corresponding to `raw_depths`.
    pub raw_depth_frequencies: Vec<u64>,
    /// Normalized observed per-bin signal before the model transformations.
    pub observed_bias: Vec<f64>,
    /// Smoothed curve from the first mixture fit.
    pub first_fitted_bias: Vec<f64>,
    /// Observed curve after noise and skew division.
    pub first_corrected_bias: Vec<f64>,
    /// Smoothed curve from the second mixture fit.
    pub second_fitted_bias: Vec<f64>,
    /// Symmetric 166 bp target curve.
    pub target_bias: Vec<f64>,
    /// Per-bin observed-to-first-fit noise division factors.
    pub noise_division_factors: Vec<f64>,
    /// Per-bin normalized linear skew division factors.
    pub skew_division_factors: Vec<f64>,
    /// Per-bin intermediate-to-target mean-shift division factors.
    pub mean_shift_division_factors: Vec<f64>,
    /// Per-bin application lookup equal to `1 / (noise * skew * mean)`.
    pub combined_weights: Vec<f64>,
    /// Parameters optimized against the original observed curve.
    pub initial_fit: MixtureFitParameters,
    /// Parameters optimized after noise and skew correction.
    pub refit: MixtureFitParameters,
    /// Minimum `forward.pos` to `reverse.reference_end` fragment length used during fitting.
    pub minimum_fragment_length: u32,
    /// Inclusive maximum `forward.pos` to `reverse.reference_end` fragment length used during
    /// fitting.
    pub maximum_fragment_length: u32,
    /// Minimum read mapping quality used during fitting.
    pub minimum_mapq: u8,
    /// Whether fitting required the BAM proper-pair flag.
    pub require_proper_pair: bool,
    /// Whether each accepted read was treated as a complete fragment.
    pub reads_are_fragments: bool,
    /// Whether reference gaps were excluded from covered positions.
    pub ignore_gap: bool,
    /// Whether blacklisted positions were excluded during fitting.
    pub blacklist_used: bool,
}

impl OverlappingLengthsCorrectionPackage {
    /// Assemble the persisted package from command settings, statistics, and fitted curves.
    ///
    /// Only scientific settings and derived values are copied. Input paths and file fingerprints
    /// are intentionally absent so a valid model can be applied to a subset BAM without retaining
    /// sensitive source locations.
    pub(crate) fn from_model(
        config: &OverlappingLengthsCorrectionConfig,
        statistics: &OverlappingLengthStatistics,
        model: OverlappingLengthModel,
    ) -> Self {
        let (raw_depths, raw_depth_frequencies): (Vec<_>, Vec<_>) = statistics
            .raw_depth_frequencies
            .iter()
            .map(|(&depth, &frequency)| (depth, frequency))
            .unzip();
        Self {
            version: OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION,
            length_bin_edges: statistics.length_bin_edges.clone(),
            length_bin_midpoints: model.bin_midpoints,
            length_bin_base_counts: statistics.length_bin_base_counts.clone(),
            observed_signal_sums: statistics.observed_signal_sums.clone(),
            raw_depths,
            raw_depth_frequencies,
            observed_bias: model.observed_bias,
            first_fitted_bias: model.first_fitted_bias,
            first_corrected_bias: model.first_corrected_bias,
            second_fitted_bias: model.second_fitted_bias,
            target_bias: model.target_bias,
            noise_division_factors: model.noise_division_factors,
            skew_division_factors: model.skew_division_factors,
            mean_shift_division_factors: model.mean_shift_division_factors,
            combined_weights: model.combined_weights,
            initial_fit: model.initial_fit,
            refit: model.refit,
            minimum_fragment_length: config.min_fragment_length,
            maximum_fragment_length: config.max_fragment_length,
            minimum_mapq: config.min_mapq,
            require_proper_pair: config.require_proper_pair,
            reads_are_fragments: config.unpaired.reads_are_fragments,
            ignore_gap: config.ignore_gap,
            blacklist_used: config.blacklist.is_some(),
        }
    }

    /// Write the complete model package as a Zarr directory.
    ///
    /// Every curve uses the `length_bin` dimension. The raw-depth values and frequencies share the
    /// `raw_depth` dimension, while bin edges use `length_edge`. Root attributes store algorithm
    /// constants and compatibility settings but never source paths.
    pub fn write_zarr(&self, path: impl AsRef<Path>) -> Result<()> {
        self.validate()?;
        let store = create_zarr_store(path.as_ref(), "overlapping fragment length model")?;
        write_zarr_root_metadata(
            store.clone(),
            "overlapping fragment length model",
            json!({
                "cfdnalab_schema": OVERLAPPING_LENGTHS_CORRECTION_SCHEMA,
                "cfdnalab_schema_version": self.version,
                "package_role": "sample_overlapping_fragment_length_normalization_model",
                "correction_units": "multiplicative_positional_coverage_weight",
                "fragment_length_definition": "forward.pos_to_reverse.reference_end",
                "validated_fragment_length_range": [100, 220],
                "base_sigma": BASE_SIGMA,
                "student_t_degrees_of_freedom": STUDENT_T_DEGREES_OF_FREEDOM,
                "target_mean_fragment_length": TARGET_MEAN_FRAGMENT_LENGTH,
                "minimum_fragment_length": self.minimum_fragment_length,
                "maximum_fragment_length": self.maximum_fragment_length,
                "minimum_mapq": self.minimum_mapq,
                "require_proper_pair": self.require_proper_pair,
                "reads_are_fragments": self.reads_are_fragments,
                "ignore_gap": self.ignore_gap,
                "blacklist_used": self.blacklist_used,
                "initial_fit": fit_json(self.initial_fit),
                "refit": fit_json(self.refit),
            }),
        )?;

        write_array_f64(
            store.clone(),
            "length_bin_edges",
            "length_edge",
            &self.length_bin_edges,
        )?;
        write_array_f64(
            store.clone(),
            "length_bin_midpoints",
            "length_bin",
            &self.length_bin_midpoints,
        )?;
        write_array_u64(
            store.clone(),
            "length_bin_base_counts",
            "length_bin",
            &self.length_bin_base_counts,
        )?;
        write_array_f64(
            store.clone(),
            "observed_signal_sums",
            "length_bin",
            &self.observed_signal_sums,
        )?;
        write_array_u32(store.clone(), "raw_depths", "raw_depth", &self.raw_depths)?;
        write_array_u64(
            store.clone(),
            "raw_depth_frequencies",
            "raw_depth",
            &self.raw_depth_frequencies,
        )?;
        // All fitted and normalization curves have exactly one value per average-length bin
        for (name, values) in [
            ("observed_bias", &self.observed_bias),
            ("first_fitted_bias", &self.first_fitted_bias),
            ("first_corrected_bias", &self.first_corrected_bias),
            ("second_fitted_bias", &self.second_fitted_bias),
            ("target_bias", &self.target_bias),
            ("noise_division_factors", &self.noise_division_factors),
            ("skew_division_factors", &self.skew_division_factors),
            (
                "mean_shift_division_factors",
                &self.mean_shift_division_factors,
            ),
            ("combined_weights", &self.combined_weights),
        ] {
            write_array_f64(store.clone(), name, "length_bin", values)?;
        }
        Ok(())
    }

    /// Load and validate an overlapping fragment length model Zarr package.
    ///
    /// Schema identity, schema version, required attributes, array lengths, and application weights
    /// are checked before the package can influence fcoverage.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        if !path.is_dir() {
            bail!(
                "overlapping fragment length model must be an existing .zarr directory: {}",
                path.display()
            );
        }
        let root = read_zarr_root_attributes(path)?;
        ensure_zarr_schema(
            &root,
            OVERLAPPING_LENGTHS_CORRECTION_SCHEMA,
            OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION,
            "overlapping fragment length model",
        )?;
        let store = Arc::new(FilesystemStore::new(path)?);
        let package = Self {
            version: u32_attribute(&root, "cfdnalab_schema_version")?,
            length_bin_edges: read_zarr_array1(store.clone(), "/length_bin_edges")?,
            length_bin_midpoints: read_zarr_array1(store.clone(), "/length_bin_midpoints")?,
            length_bin_base_counts: read_zarr_array1(store.clone(), "/length_bin_base_counts")?,
            observed_signal_sums: read_zarr_array1(store.clone(), "/observed_signal_sums")?,
            raw_depths: read_zarr_array1(store.clone(), "/raw_depths")?,
            raw_depth_frequencies: read_zarr_array1(store.clone(), "/raw_depth_frequencies")?,
            observed_bias: read_zarr_array1(store.clone(), "/observed_bias")?,
            first_fitted_bias: read_zarr_array1(store.clone(), "/first_fitted_bias")?,
            first_corrected_bias: read_zarr_array1(store.clone(), "/first_corrected_bias")?,
            second_fitted_bias: read_zarr_array1(store.clone(), "/second_fitted_bias")?,
            target_bias: read_zarr_array1(store.clone(), "/target_bias")?,
            noise_division_factors: read_zarr_array1(store.clone(), "/noise_division_factors")?,
            skew_division_factors: read_zarr_array1(store.clone(), "/skew_division_factors")?,
            mean_shift_division_factors: read_zarr_array1(
                store.clone(),
                "/mean_shift_division_factors",
            )?,
            combined_weights: read_zarr_array1(store, "/combined_weights")?,
            initial_fit: fit_from_json(&root, "initial_fit")?,
            refit: fit_from_json(&root, "refit")?,
            minimum_fragment_length: u32_attribute(&root, "minimum_fragment_length")?,
            maximum_fragment_length: u32_attribute(&root, "maximum_fragment_length")?,
            minimum_mapq: u8::try_from(u32_attribute(&root, "minimum_mapq")?)
                .context("minimum_mapq must fit in u8")?,
            require_proper_pair: bool_attribute(&root, "require_proper_pair")?,
            reads_are_fragments: bool_attribute(&root, "reads_are_fragments")?,
            ignore_gap: bool_attribute(&root, "ignore_gap")?,
            blacklist_used: bool_attribute(&root, "blacklist_used")?,
        };
        package.validate()?;
        Ok(package)
    }

    /// Map a positional average fragment length to its clipped package bin index.
    ///
    /// Returning the index separately lets tiled inference retain compact integer values and defer
    /// multiplier lookup until finalized coverage is available. Values outside the fitted interval
    /// use the nearest extreme bin, matching LIONHEART's clipped `numpy.digitize` lookup.
    pub(crate) fn bin_index_for_average_length(&self, average_length: f64) -> Result<usize> {
        ensure!(
            average_length.is_finite() && average_length > 0.0,
            "average overlapping fragment length must be finite and positive"
        );
        let insertion = self
            .length_bin_edges
            .partition_point(|edge| *edge <= average_length);
        let bin_index = insertion
            .saturating_sub(1)
            .min(self.combined_weights.len() - 1);
        Ok(bin_index)
    }

    /// Read a combined positional multiplier by its persisted length-bin index.
    pub(crate) fn weight_for_bin_index(&self, bin_index: u32) -> Result<f64> {
        self.combined_weights
            .get(bin_index as usize)
            .copied()
            .with_context(|| {
                format!(
                    "overlap-length bin index {} is outside the {} package bins",
                    bin_index,
                    self.combined_weights.len()
                )
            })
    }

    /// Validate structural invariants required by fitting output and lookup application.
    fn validate(&self) -> Result<()> {
        ensure!(
            self.version == OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION,
            "unsupported overlapping fragment length model schema version {}",
            self.version
        );
        ensure!(
            self.length_bin_edges.len() >= 2
                && self
                    .length_bin_edges
                    .windows(2)
                    .all(|pair| pair[0] < pair[1]),
            "overlapping fragment length bin edges must increase strictly"
        );
        let num_bins = self.length_bin_edges.len() - 1;
        for (name, length) in [
            ("length_bin_midpoints", self.length_bin_midpoints.len()),
            ("length_bin_base_counts", self.length_bin_base_counts.len()),
            ("observed_signal_sums", self.observed_signal_sums.len()),
            ("observed_bias", self.observed_bias.len()),
            ("first_fitted_bias", self.first_fitted_bias.len()),
            ("first_corrected_bias", self.first_corrected_bias.len()),
            ("second_fitted_bias", self.second_fitted_bias.len()),
            ("target_bias", self.target_bias.len()),
            ("noise_division_factors", self.noise_division_factors.len()),
            ("skew_division_factors", self.skew_division_factors.len()),
            (
                "mean_shift_division_factors",
                self.mean_shift_division_factors.len(),
            ),
            ("combined_weights", self.combined_weights.len()),
        ] {
            ensure!(
                length == num_bins,
                "{} has {} values but expected {}",
                name,
                length,
                num_bins
            );
        }
        ensure!(
            self.raw_depths.len() == self.raw_depth_frequencies.len(),
            "raw depth values and frequencies have different lengths"
        );
        ensure!(
            self.combined_weights
                .iter()
                .all(|weight| weight.is_finite() && *weight > 0.0),
            "combined overlapping fragment length weights must be finite and positive"
        );
        ensure!(
            self.minimum_fragment_length <= self.maximum_fragment_length,
            "package fragment length interval is invalid"
        );
        Ok(())
    }
}

/// Convert fitted parameters to a stable root-attribute object.
fn fit_json(parameters: MixtureFitParameters) -> Value {
    json!({
        "scale_multiplier": parameters.scale_multiplier,
        "skewness": parameters.skewness,
        "mean_fragment_length": parameters.mean_fragment_length,
    })
}

/// Parse a named fitted-parameter object from Zarr root attributes.
fn fit_from_json(root: &Value, key: &str) -> Result<MixtureFitParameters> {
    let fit = root
        .get(key)
        .and_then(Value::as_object)
        .with_context(|| format!("missing {key} metadata"))?;
    let number = |field: &str| {
        fit.get(field)
            .and_then(Value::as_f64)
            .with_context(|| format!("missing {key}.{field} metadata"))
    };
    Ok(MixtureFitParameters {
        scale_multiplier: number("scale_multiplier")?,
        skewness: number("skewness")?,
        mean_fragment_length: number("mean_fragment_length")?,
    })
}

/// Read a required non-negative integer root attribute as `u32`.
fn u32_attribute(root: &Value, key: &str) -> Result<u32> {
    u32::try_from(
        root.get(key)
            .and_then(Value::as_u64)
            .with_context(|| format!("missing {key} metadata"))?,
    )
    .with_context(|| format!("{key} metadata must fit in u32"))
}

/// Read a required Boolean root attribute.
fn bool_attribute(root: &Value, key: &str) -> Result<bool> {
    root.get(key)
        .and_then(Value::as_bool)
        .with_context(|| format!("missing {key} metadata"))
}

/// Write a single-chunk floating-point array with one named dimension.
fn write_array_f64(
    store: Arc<FilesystemStore>,
    name: &str,
    dimension: &str,
    values: &[f64],
) -> Result<()> {
    write_single_chunk_zarr_array(
        store,
        name,
        &[values.len()],
        &[dimension],
        values,
        data_type::float64(),
        0.0,
        json!({"long_name": name}),
    )
}

/// Write a single-chunk `u32` array with one named dimension.
fn write_array_u32(
    store: Arc<FilesystemStore>,
    name: &str,
    dimension: &str,
    values: &[u32],
) -> Result<()> {
    write_single_chunk_zarr_array(
        store,
        name,
        &[values.len()],
        &[dimension],
        values,
        data_type::uint32(),
        0u32,
        json!({"long_name": name}),
    )
}

/// Write a single-chunk `u64` array with one named dimension.
fn write_array_u64(
    store: Arc<FilesystemStore>,
    name: &str,
    dimension: &str,
    values: &[u64],
) -> Result<()> {
    write_single_chunk_zarr_array(
        store,
        name,
        &[values.len()],
        &[dimension],
        values,
        data_type::uint64(),
        0u64,
        json!({"long_name": name}),
    )
}

#[cfg(test)]
mod tests {
    include!("package_tests.rs");
}
