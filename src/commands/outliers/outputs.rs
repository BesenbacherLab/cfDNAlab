//! Text and plot outputs for `cfdna outliers`.
//!
//! The command module owns BAM processing and region calling. This module keeps output paths,
//! schemas, metadata, and final-file registration together so the orchestration remains readable.

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use fxhash::FxHashMap;

use crate::{
    commands::outliers::{
        config::OutliersConfig,
        distribution::{TwoStageCoverageModel, diagnose_fit},
        model::CoverageCounts,
        striding::{CoreHistogram, CoreModel},
    },
    shared::{
        bam::Contigs,
        io::{FinalOutputFiles, dot_join},
        reference::ContigFootprintEntry,
        writers::open_zstd_auto_writer,
    },
};

const HISTOGRAM_ZSTD_LEVEL: i32 = 3;

#[cfg(feature = "plotters")]
use super::plotting::write_global_fit_plot;

#[derive(Debug)]
pub(super) struct OutputPaths {
    pub(super) keep_weights: PathBuf,
    pub(super) exact_blacklist: PathBuf,
    pub(super) flanked_blacklist: PathBuf,
    pub(super) histograms: PathBuf,
    pub(super) models: PathBuf,
    pub(super) global_fit: PathBuf,
    #[cfg(feature = "plotters")]
    pub(super) plot: PathBuf,
    pub(super) audit: PathBuf,
    pub(super) summary: PathBuf,
    pub(super) preserved_regions: PathBuf,
    pub(super) readme: PathBuf,
}

impl OutputPaths {
    pub(super) fn new(output_dir: &Path, prefix: &str) -> Self {
        Self {
            keep_weights: output_dir.join(dot_join(&[prefix, "outliers.keep_weights.tsv"])),
            exact_blacklist: output_dir.join(dot_join(&[prefix, "outliers.exact.bed"])),
            flanked_blacklist: output_dir.join(dot_join(&[prefix, "outliers.flanked.bed"])),
            histograms: output_dir.join(dot_join(&[prefix, "outliers.histograms.tsv.zst"])),
            models: output_dir.join(dot_join(&[prefix, "outliers.models.tsv"])),
            global_fit: output_dir.join(dot_join(&[prefix, "outliers.global_fit.tsv"])),
            #[cfg(feature = "plotters")]
            plot: output_dir.join(dot_join(&[prefix, "outliers.global_fit.png"])),
            audit: output_dir.join(dot_join(&[prefix, "outliers.audit.tsv"])),
            summary: output_dir.join(dot_join(&[prefix, "outliers.summary.tsv"])),
            preserved_regions: output_dir
                .join(dot_join(&[prefix, "outliers.preserved_regions.tsv"])),
            readme: output_dir.join(dot_join(&[prefix, "outliers.README.txt"])),
        }
    }

    pub(super) fn all(&self) -> Vec<PathBuf> {
        let mut paths = vec![
            self.keep_weights.clone(),
            self.exact_blacklist.clone(),
            self.flanked_blacklist.clone(),
            self.histograms.clone(),
            self.models.clone(),
            self.global_fit.clone(),
        ];
        #[cfg(feature = "plotters")]
        paths.push(self.plot.clone());
        paths.extend([
            self.audit.clone(),
            self.summary.clone(),
            self.preserved_regions.clone(),
            self.readme.clone(),
        ]);
        paths
    }

    /// Give each optional preserved-region PNG the same prefix as its summary table.
    pub(super) fn preserved_plot(&self, region_id: &str) -> PathBuf {
        self.preserved_regions
            .with_extension(format!("{region_id}.png"))
    }
}

/// Scientifically relevant settings shared by every text output.
pub(super) struct OutputMetadata<'a> {
    pub(super) config: &'a OutliersConfig,
    pub(super) chromosomes: &'a [String],
    pub(super) contigs: &'a Contigs,
    pub(super) tail_probability: f64,
    pub(super) eligible_positions: u64,
    pub(super) blacklist_flank: u32,
}

pub(super) fn write_common_metadata(
    writer: &mut impl Write,
    metadata: &OutputMetadata<'_>,
) -> Result<()> {
    let config = metadata.config;
    writeln!(writer, "# cfdnalab_version={}", env!("CARGO_PKG_VERSION"))?;
    writeln!(writer, "# command=cfdna outliers")?;
    writeln!(
        writer,
        "# input_alignment={}",
        serde_json::to_string(&config.ioc.bam.to_string_lossy())?
    )?;
    writeln!(writer, "# coordinate_system=0-based-half-open")?;
    writeln!(writer, "# detector=two_stage_{}", config.model.as_str())?;
    writeln!(writer, "# model={}", config.model.as_str())?;
    writeln!(writer, "# local_mean_target=active_component_mean")?;
    writeln!(
        writer,
        "# protection_window_size={}",
        config.protection_window_size
    )?;
    writeln!(
        writer,
        "# min_protection_positions={}",
        config.min_protection_positions
    )?;
    writeln!(
        writer,
        "# max_outlier_percentage={}",
        config.max_outlier_percentage
    )?;
    writeln!(
        writer,
        "# ambiguous_region_action={}",
        config.ambiguous_region_action.as_str()
    )?;
    let mut footprint: Vec<_> = metadata
        .contigs
        .contigs
        .iter()
        .map(|(name, &(_, length))| ContigFootprintEntry {
            name: name.clone(),
            size: length as u64,
        })
        .collect();
    footprint.sort_unstable_by(|left, right| left.name.cmp(&right.name));
    writeln!(
        writer,
        "# analysed_chromosomes={}",
        serde_json::to_string(&footprint)?
    )?;
    writeln!(writer, "# target={}", config.target.as_str())?;
    writeln!(
        writer,
        "# plotting_available={}",
        cfg!(feature = "plotters")
    )?;
    writeln!(
        writer,
        "# plot_preserved_regions={}",
        config.plot_preserved_regions
    )?;
    writeln!(
        writer,
        "# preserved_plot_context_multiplier={}",
        config.preserved_plot_context_multiplier
    )?;
    writeln!(
        writer,
        "# preserved_plot_min_span={}",
        config.preserved_plot_min_span
    )?;
    writeln!(writer, "# tail_probability={}", metadata.tail_probability)?;
    writeln!(
        writer,
        "# tail_probability_mode={}",
        if config.tail_probability.is_some() {
            "manual"
        } else {
            "automatic"
        }
    )?;
    writeln!(
        writer,
        "# tail_probability_multiplier={}",
        if config.tail_probability.is_none() {
            config.tail_probability_multiplier.to_string()
        } else {
            "NA".to_string()
        }
    )?;
    writeln!(
        writer,
        "# automatic_probability_denominator={}",
        if config.tail_probability.is_none() {
            metadata.eligible_positions.to_string()
        } else {
            "NA".to_string()
        }
    )?;
    writeln!(
        writer,
        "# eligible_positions={}",
        metadata.eligible_positions
    )?;
    writeln!(writer, "# model_core_size={}", config.stride)?;
    writeln!(writer, "# minimum_context_span={}", config.bin_size)?;
    writeln!(
        writer,
        "# minimum_context_fragment_support={}",
        config.min_context_fragments
    )?;
    writeln!(writer, "# tile_size={}", config.tile_size)?;
    writeln!(
        writer,
        "# raw_coverage_source=uncorrected_mapped_reference_segments"
    )?;
    writeln!(writer, "# ignore_gap={}", config.ignore_gap)?;
    writeln!(
        writer,
        "# reads_are_fragments={}",
        config.unpaired.reads_are_fragments
    )?;
    writeln!(
        writer,
        "# require_proper_pair={}",
        config.require_proper_pair
    )?;
    writeln!(writer, "# minimum_mapping_quality={}", config.min_mapq)?;
    writeln!(
        writer,
        "# fragment_length_range={}-{}",
        config.fragment_lengths.min_fragment_length, config.fragment_lengths.max_fragment_length
    )?;
    writeln!(writer, "# blacklist_flank={}", metadata.blacklist_flank)?;
    let blacklist_paths = config.blacklist.as_deref().unwrap_or_default();
    writeln!(
        writer,
        "# blacklists={}",
        serde_json::to_string(
            &blacklist_paths
                .iter()
                .map(|path| path.to_string_lossy())
                .collect::<Vec<_>>()
        )?
    )?;
    Ok(())
}

/// Write per-core models and the complete global histogram, plus a plot when plotting is compiled.
pub(super) fn write_diagnostic_outputs(
    metadata: &OutputMetadata<'_>,
    global_histogram: &CoreHistogram,
    global_model: Option<TwoStageCoverageModel>,
    models_by_chromosome: &FxHashMap<String, Vec<Option<CoreModel>>>,
    output_paths: &OutputPaths,
    final_outputs: &mut FinalOutputFiles,
) -> Result<()> {
    let models_temp = final_outputs.temp_path_for(&output_paths.models)?;
    let mut models_writer =
        BufWriter::new(File::create(&models_temp).context("creating outlier model-summary TSV")?);
    write_common_metadata(&mut models_writer, metadata)?;
    writeln!(
        models_writer,
        "chromosome\tstart\tend\tcontext_start\tcontext_end\tmodel_source\teligible_positions\tfragments_by_midpoint\tinitial_component_mean\tinitial_zero_inflation\tinitial_threshold\tsecond_fit_retained_positions\tunderlying_component_mean\tunderlying_zero_inflation\tunderlying_mean\tfinal_threshold\tobserved_positive_coverage_mean\tobserved_positive_coverage_variance\tunderlying_expected_positive_coverage_variance\tpositive_coverage_variance_ratio\tobserved_initial_tail_positions\texpected_initial_tail_positions\tobserved_final_tail_positions\texpected_final_tail_positions\tinitial_dispersion\tunderlying_dispersion\tinitial_poisson_boundary\tunderlying_poisson_boundary\tfinal_threshold_exceeds_initial\tinitial_optimizer_diagnostics\tunderlying_optimizer_diagnostics"
    )?;
    for chromosome in metadata.chromosomes {
        let models = models_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing models for chromosome '{chromosome}'"))?;
        for model in models.iter().flatten() {
            write!(
                models_writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chromosome,
                model.core.start(),
                model.core.end(),
                model.context.start(),
                model.context.end(),
                model.source.as_str(),
                model.eligible_positions,
                model.fragment_support,
                model.model.initial.component_mean,
                model.model.initial.zero_inflation,
                model.model.initial_threshold,
                model.model.second_fit_retained_positions,
                model.model.underlying_fit.component_mean,
                model.model.underlying_fit.zero_inflation,
                model.model.underlying_fit.mean(),
                model.model.final_threshold,
                model.diagnostics.observed_positive_coverage_mean,
                model.diagnostics.observed_positive_coverage_variance,
                model
                    .diagnostics
                    .underlying_expected_positive_coverage_variance,
                model.diagnostics.positive_coverage_variance_ratio,
                model.diagnostics.observed_initial_tail_positions,
                model.diagnostics.expected_initial_tail_positions,
                model.diagnostics.observed_final_tail_positions,
                model.diagnostics.expected_final_tail_positions,
                model.model.initial.dispersion,
                model.model.underlying_fit.dispersion,
                model.model.initial.dispersion == 0.0,
                model.model.underlying_fit.dispersion == 0.0,
                model.model.final_threshold > model.model.initial_threshold
            )?;
            writeln!(
                models_writer,
                "\t{}\t{}",
                model
                    .model
                    .initial_optimizer
                    .map(|diagnostics| serde_json::to_string(&diagnostics))
                    .transpose()?
                    .unwrap_or_else(|| "NA".to_string()),
                model
                    .model
                    .underlying_optimizer
                    .map(|diagnostics| serde_json::to_string(&diagnostics))
                    .transpose()?
                    .unwrap_or_else(|| "NA".to_string())
            )?;
        }
    }
    models_writer.flush()?;
    drop(models_writer);
    final_outputs.record(models_temp, output_paths.models.clone())?;

    let global_fit_temp = final_outputs.temp_path_for(&output_paths.global_fit)?;
    let mut global_writer =
        BufWriter::new(File::create(&global_fit_temp).context("creating global outlier fit TSV")?);
    write_common_metadata(&mut global_writer, metadata)?;
    writeln!(
        global_writer,
        "# global_fit_status={}",
        if global_model.is_some() {
            "required_for_fallback"
        } else {
            "not_needed"
        }
    )?;
    if let Some(global_model) = global_model {
        writeln!(
            global_writer,
            "# initial_component_mean={}",
            global_model.initial.component_mean
        )?;
        writeln!(
            global_writer,
            "# initial_zero_inflation={}",
            global_model.initial.zero_inflation
        )?;
        writeln!(
            global_writer,
            "# initial_threshold={}",
            global_model.initial_threshold
        )?;
        writeln!(
            global_writer,
            "# underlying_component_mean={}",
            global_model.underlying_fit.component_mean
        )?;
        writeln!(
            global_writer,
            "# underlying_zero_inflation={}",
            global_model.underlying_fit.zero_inflation
        )?;
        writeln!(
            global_writer,
            "# underlying_mean={}",
            global_model.underlying_fit.mean()
        )?;
        writeln!(
            global_writer,
            "# final_threshold={}",
            global_model.final_threshold
        )?;
        writeln!(
            global_writer,
            "# initial_dispersion={}",
            global_model.initial.dispersion
        )?;
        writeln!(
            global_writer,
            "# underlying_dispersion={}",
            global_model.underlying_fit.dispersion
        )?;
        writeln!(
            global_writer,
            "# initial_poisson_boundary={}",
            global_model.initial.dispersion == 0.0
        )?;
        writeln!(
            global_writer,
            "# underlying_poisson_boundary={}",
            global_model.underlying_fit.dispersion == 0.0
        )?;
        writeln!(
            global_writer,
            "# initial_optimizer_diagnostics={}",
            serde_json::to_string(&global_model.initial_optimizer)?
        )?;
        writeln!(
            global_writer,
            "# underlying_optimizer_diagnostics={}",
            serde_json::to_string(&global_model.underlying_optimizer)?
        )?;
        let global_diagnostics = diagnose_fit(&global_histogram.counts, global_model)?;
        writeln!(
            global_writer,
            "# observed_positive_coverage_mean={}",
            global_diagnostics.observed_positive_coverage_mean
        )?;
        writeln!(
            global_writer,
            "# observed_positive_coverage_variance={}",
            global_diagnostics.observed_positive_coverage_variance
        )?;
        writeln!(
            global_writer,
            "# underlying_expected_positive_coverage_variance={}",
            global_diagnostics.underlying_expected_positive_coverage_variance
        )?;
        writeln!(
            global_writer,
            "# positive_coverage_variance_ratio={}",
            global_diagnostics.positive_coverage_variance_ratio
        )?;
        writeln!(
            global_writer,
            "# observed_initial_tail_positions={}",
            global_diagnostics.observed_initial_tail_positions
        )?;
        writeln!(
            global_writer,
            "# expected_initial_tail_positions={}",
            global_diagnostics.expected_initial_tail_positions
        )?;
        writeln!(
            global_writer,
            "# observed_final_tail_positions={}",
            global_diagnostics.observed_final_tail_positions
        )?;
        writeln!(
            global_writer,
            "# expected_final_tail_positions={}",
            global_diagnostics.expected_final_tail_positions
        )?;
    } else {
        writeln!(
            global_writer,
            "# unavailable_fit_columns=NA because no core needs a global fallback model"
        )?;
    }
    writeln!(
        global_writer,
        "coverage\tobserved_positions\tinitial_fitted_positions\tunderlying_fitted_positions"
    )?;
    let global_positions = global_histogram.eligible_positions()? as f64;
    let maximum_observed_coverage = maximum_observed_coverage(&global_histogram.counts)?;
    for coverage in 0..=maximum_observed_coverage {
        let observed = global_histogram.counts.get(&coverage).copied().unwrap_or(0);
        if let Some(model) = global_model {
            writeln!(
                global_writer,
                "{coverage}\t{observed}\t{}\t{}",
                global_positions * model.initial.probability_mass(coverage),
                global_positions * model.underlying_fit.probability_mass(coverage)
            )?;
        } else {
            writeln!(global_writer, "{coverage}\t{observed}\tNA\tNA")?;
        }
    }
    global_writer.flush()?;
    drop(global_writer);
    final_outputs.record(global_fit_temp, output_paths.global_fit.clone())?;

    #[cfg(feature = "plotters")]
    {
        let plot_temp = final_outputs.temp_path_for(&output_paths.plot)?;
        write_global_fit_plot(
            &plot_temp,
            &global_histogram.counts,
            global_model,
            metadata.chromosomes,
            models_by_chromosome,
        )?;
        final_outputs.record(plot_temp, output_paths.plot.clone())?;
    }
    Ok(())
}

/// Write sparse per-core coverage histograms with midpoint-assigned fragment support.
pub(super) fn write_histogram_output(
    metadata: &OutputMetadata<'_>,
    histograms_by_chromosome: &FxHashMap<String, Vec<CoreHistogram>>,
    output_path: &Path,
    final_outputs: &mut FinalOutputFiles,
) -> Result<()> {
    let histogram_temp = final_outputs.temp_path_for(output_path)?;
    let mut writer = open_zstd_auto_writer(
        &histogram_temp,
        HISTOGRAM_ZSTD_LEVEL,
        Some(metadata.config.ioc.n_threads as u32),
    )?;
    write_common_metadata(&mut writer, metadata)?;
    writeln!(
        writer,
        "chromosome\tstart\tend\tfragments_by_midpoint\tcoverage\tobserved_positions"
    )?;
    for chromosome in metadata.chromosomes {
        let histograms = histograms_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing histograms for chromosome '{chromosome}'"))?;
        for histogram in histograms {
            if histogram.counts.is_empty() {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t0\t0",
                    chromosome,
                    histogram.interval.start(),
                    histogram.interval.end(),
                    histogram.fragment_support,
                )?;
                continue;
            }
            for (&coverage, &positions) in &histogram.counts {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    chromosome,
                    histogram.interval.start(),
                    histogram.interval.end(),
                    histogram.fragment_support,
                    coverage,
                    positions
                )?;
            }
        }
    }
    writer.flush()?;
    drop(writer);
    final_outputs.record(histogram_temp, output_path.to_path_buf())?;
    Ok(())
}

fn maximum_observed_coverage(observed_counts: &CoverageCounts) -> Result<u32> {
    observed_counts
        .last_key_value()
        .map(|(&coverage, _)| coverage)
        .context("global histogram has no observed coverage values")
}

#[cfg(test)]
mod tests {
    include!("outputs_tests.rs");
}
