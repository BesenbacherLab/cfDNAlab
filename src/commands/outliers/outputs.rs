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
        model::{CoverageCounts, TwoStageZipModel, diagnose_zip_fit},
        outliers::CalledRegion,
        plotting::write_global_fit_plot,
        striding::{CoreHistogram, CoreModel},
    },
    shared::{
        bam::Contigs,
        interval::{Interval, TouchingMergePolicy, push_merged_interval},
        io::{FinalOutputFiles, dot_join},
        writers::open_zstd_auto_writer,
    },
};

const HISTOGRAM_ZSTD_LEVEL: i32 = 3;

#[derive(Debug)]
pub(super) struct OutputPaths {
    pub(super) keep_weights: PathBuf,
    pub(super) exact_blacklist: PathBuf,
    pub(super) flanked_blacklist: PathBuf,
    pub(super) histograms: PathBuf,
    pub(super) models: PathBuf,
    pub(super) global_fit: PathBuf,
    pub(super) plot: PathBuf,
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
            plot: output_dir.join(dot_join(&[prefix, "outliers.global_fit.png"])),
        }
    }

    pub(super) fn all(&self) -> Vec<PathBuf> {
        vec![
            self.keep_weights.clone(),
            self.exact_blacklist.clone(),
            self.flanked_blacklist.clone(),
            self.histograms.clone(),
            self.models.clone(),
            self.global_fit.clone(),
            self.plot.clone(),
        ]
    }
}

/// Scientifically relevant settings shared by every text output.
pub(super) struct OutputMetadata<'a> {
    pub(super) config: &'a OutliersConfig,
    pub(super) chromosomes: &'a [String],
    pub(super) tail_probability: f64,
    pub(super) eligible_positions: u64,
    pub(super) blacklist_flank: u32,
}

fn write_common_metadata(writer: &mut impl Write, metadata: &OutputMetadata<'_>) -> Result<()> {
    let config = metadata.config;
    writeln!(writer, "# cfdnalab_version={}", env!("CARGO_PKG_VERSION"))?;
    writeln!(writer, "# command=cfdna outliers")?;
    writeln!(
        writer,
        "# input_alignment={}",
        serde_json::to_string(&config.ioc.bam.to_string_lossy())?
    )?;
    writeln!(
        writer,
        "# selected_chromosomes={}",
        serde_json::to_string(metadata.chromosomes)?
    )?;
    writeln!(writer, "# coordinate_system=0-based-half-open")?;
    writeln!(writer, "# detector=two_stage_zero_inflated_poisson")?;
    writeln!(writer, "# hysteresis=false")?;
    writeln!(writer, "# target={}", config.target.as_str())?;
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
        config.tail_probability_multiplier
    )?;
    writeln!(
        writer,
        "# automatic_probability_denominator={}",
        metadata.eligible_positions
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
        "# fragment_support_definition=accepted_fragments_with_unblacklisted_midpoint"
    )?;
    writeln!(
        writer,
        "# raw_coverage_source=uncorrected_mapped_reference_segments"
    )?;
    writeln!(writer, "# ignore_inter_mate_gap={}", config.ignore_gap)?;
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
        "# input_blacklists={}",
        serde_json::to_string(
            &blacklist_paths
                .iter()
                .map(|path| path.to_string_lossy())
                .collect::<Vec<_>>()
        )?
    )?;
    Ok(())
}

/// Explain the fit diagnostics in each output that carries them.
fn write_fit_diagnostic_definitions(writer: &mut impl Write) -> Result<()> {
    writeln!(
        writer,
        "# positive_coverage_diagnostic_scope=observed values use the complete original positive-coverage histogram, including the extreme tail, while the ZIP expectation uses the underlying model estimated by the second fit"
    )?;
    writeln!(
        writer,
        "# observed_positive_coverage_variance_definition=population mean squared deviation among original positions with coverage greater than zero"
    )?;
    writeln!(
        writer,
        "# underlying_zip_expected_positive_coverage_variance_definition=variance among positive values expected from the fitted underlying ZIP"
    )?;
    writeln!(
        writer,
        "# positive_coverage_variance_ratio_definition=observed_positive_coverage_variance divided by underlying_zip_expected_positive_coverage_variance"
    )?;
    writeln!(
        writer,
        "# positive_coverage_variance_ratio_interpretation=near 1 is Poisson-like spread, above 1 is more variable, and below 1 is less variable. This diagnostic does not affect fitting or calling"
    )?;
    writeln!(
        writer,
        "# tail_diagnostic_definition=observed and expected positions at or above T1 and T2. These are more direct checks of calling-tail calibration"
    )?;
    writeln!(
        writer,
        "# unavailable_positive_coverage_diagnostics=NaN when an all-zero incomplete context uses the global fallback model"
    )?;
    Ok(())
}

/// Write per-core models, the complete global observed and fitted histogram, and the combined plot.
pub(super) fn write_diagnostic_outputs(
    metadata: &OutputMetadata<'_>,
    global_histogram: &CoreHistogram,
    global_model: TwoStageZipModel,
    models_by_chromosome: &FxHashMap<String, Vec<CoreModel>>,
    output_paths: &OutputPaths,
    final_outputs: &mut FinalOutputFiles,
) -> Result<()> {
    let models_temp = final_outputs.temp_path_for(&output_paths.models)?;
    let mut models_writer =
        BufWriter::new(File::create(&models_temp).context("creating outlier model-summary TSV")?);
    write_common_metadata(&mut models_writer, metadata)?;
    write_fit_diagnostic_definitions(&mut models_writer)?;
    writeln!(
        models_writer,
        "chromosome\tstart\tend\tcontext_start\tcontext_end\tmodel_source\teligible_positions\tfragment_support\tinitial_lambda\tinitial_zero_inflation\tinitial_threshold\tsecond_fit_retained_positions\tunderlying_lambda\tunderlying_zero_inflation\tunderlying_mean\tfinal_threshold\tobserved_positive_coverage_mean\tobserved_positive_coverage_variance\tunderlying_zip_expected_positive_coverage_variance\tpositive_coverage_variance_ratio\tobserved_initial_tail_positions\texpected_initial_tail_positions\tobserved_final_tail_positions\texpected_final_tail_positions"
    )?;
    for chromosome in metadata.chromosomes {
        let models = models_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing models for chromosome '{chromosome}'"))?;
        for model in models {
            writeln!(
                models_writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chromosome,
                model.core.start(),
                model.core.end(),
                model.context.start(),
                model.context.end(),
                model.source.as_str(),
                model.eligible_positions,
                model.fragment_support,
                model.zip.initial.lambda,
                model.zip.initial.zero_inflation,
                model.zip.initial_threshold,
                model.zip.second_fit_retained_positions,
                model.zip.underlying_fit.lambda,
                model.zip.underlying_fit.zero_inflation,
                model.zip.underlying_fit.mean(),
                model.zip.final_threshold,
                model.diagnostics.observed_positive_coverage_mean,
                model.diagnostics.observed_positive_coverage_variance,
                model
                    .diagnostics
                    .underlying_zip_expected_positive_coverage_variance,
                model.diagnostics.positive_coverage_variance_ratio,
                model.diagnostics.observed_initial_tail_positions,
                model.diagnostics.expected_initial_tail_positions,
                model.diagnostics.observed_final_tail_positions,
                model.diagnostics.expected_final_tail_positions
            )?;
        }
    }
    models_writer.flush()?;
    drop(models_writer);
    final_outputs.record(models_temp, output_paths.models.clone())?;

    let global_fit_temp = final_outputs.temp_path_for(&output_paths.global_fit)?;
    let mut global_writer =
        BufWriter::new(File::create(&global_fit_temp).context("creating global outlier ZIP TSV")?);
    write_common_metadata(&mut global_writer, metadata)?;
    write_fit_diagnostic_definitions(&mut global_writer)?;
    writeln!(
        global_writer,
        "# initial_lambda={}",
        global_model.initial.lambda
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
        "# underlying_lambda={}",
        global_model.underlying_fit.lambda
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
    let global_diagnostics = diagnose_zip_fit(&global_histogram.counts, global_model)?;
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
        "# underlying_zip_expected_positive_coverage_variance={}",
        global_diagnostics.underlying_zip_expected_positive_coverage_variance
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
    writeln!(
        global_writer,
        "coverage\tobserved_positions\tinitial_fitted_positions\tunderlying_fitted_positions"
    )?;
    let global_positions = global_histogram.eligible_positions()? as f64;
    let maximum_observed_coverage = maximum_observed_coverage(&global_histogram.counts)?;
    for coverage in 0..=maximum_observed_coverage {
        let observed = global_histogram.counts.get(&coverage).copied().unwrap_or(0);
        writeln!(
            global_writer,
            "{}\t{}\t{}\t{}",
            coverage,
            observed,
            global_positions * global_model.initial.probability_mass(coverage),
            global_positions * global_model.underlying_fit.probability_mass(coverage)
        )?;
    }
    global_writer.flush()?;
    drop(global_writer);
    final_outputs.record(global_fit_temp, output_paths.global_fit.clone())?;

    let plot_temp = final_outputs.temp_path_for(&output_paths.plot)?;
    write_global_fit_plot(
        &plot_temp,
        &global_histogram.counts,
        global_model,
        metadata.chromosomes,
        models_by_chromosome,
    )?;
    final_outputs.record(plot_temp, output_paths.plot.clone())?;
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
        "chromosome\tstart\tend\tfragment_support\tcoverage\tobserved_positions"
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

/// Write sparse regional weights and exact and flanked called-region BED files.
pub(super) fn write_call_outputs(
    metadata: &OutputMetadata<'_>,
    contigs: &Contigs,
    regions_by_chromosome: &FxHashMap<String, Vec<CalledRegion>>,
    output_paths: &OutputPaths,
    final_outputs: &mut FinalOutputFiles,
) -> Result<()> {
    let keep_weights_temp = final_outputs.temp_path_for(&output_paths.keep_weights)?;
    let exact_blacklist_temp = final_outputs.temp_path_for(&output_paths.exact_blacklist)?;
    let flanked_blacklist_temp = final_outputs.temp_path_for(&output_paths.flanked_blacklist)?;
    let mut keep_writer = BufWriter::new(
        File::create(&keep_weights_temp).context("creating outlier keep-weight TSV")?,
    );
    let mut exact_writer =
        BufWriter::new(File::create(&exact_blacklist_temp).context("creating exact outlier BED")?);
    let mut flanked_writer = BufWriter::new(
        File::create(&flanked_blacklist_temp).context("creating flanked outlier BED")?,
    );

    write_common_metadata(&mut keep_writer, metadata)?;
    writeln!(keep_writer, "# omitted_keep_weight=1.0")?;
    writeln!(keep_writer, "# fragment_overlap_rule=minimum_keep_weight")?;
    writeln!(keep_writer, "chromosome\tstart\tend\tkeep_weight")?;
    write_common_metadata(&mut exact_writer, metadata)?;
    writeln!(exact_writer, "# output=exact_called_regions")?;
    write_common_metadata(&mut flanked_writer, metadata)?;
    writeln!(flanked_writer, "# output=flanked_called_regions")?;

    for chromosome in metadata.chromosomes {
        let regions = regions_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing called regions for chromosome '{chromosome}'"))?;
        for region in regions {
            writeln!(
                keep_writer,
                "{}\t{}\t{}\t{}",
                chromosome,
                region.interval.start(),
                region.interval.end(),
                region.keep_weight()
            )?;
            writeln!(
                exact_writer,
                "{}\t{}\t{}",
                chromosome,
                region.interval.start(),
                region.interval.end()
            )?;
        }

        let &(_, chromosome_length) = contigs
            .contigs
            .get(chromosome)
            .with_context(|| format!("missing contig metadata for '{chromosome}'"))?;
        let mut flanked_regions = Vec::with_capacity(regions.len());
        for region in regions {
            let start = region
                .interval
                .start()
                .saturating_sub(metadata.blacklist_flank);
            let end = region
                .interval
                .end()
                .saturating_add(metadata.blacklist_flank)
                .min(chromosome_length);
            push_merged_interval(
                &mut flanked_regions,
                Interval::new(start, end)?,
                TouchingMergePolicy::MergeTouching,
            );
        }
        for interval in flanked_regions {
            writeln!(
                flanked_writer,
                "{}\t{}\t{}",
                chromosome,
                interval.start(),
                interval.end()
            )?;
        }
    }

    keep_writer.flush()?;
    exact_writer.flush()?;
    flanked_writer.flush()?;
    drop((keep_writer, exact_writer, flanked_writer));
    final_outputs.record(keep_weights_temp, output_paths.keep_weights.clone())?;
    final_outputs.record(exact_blacklist_temp, output_paths.exact_blacklist.clone())?;
    final_outputs.record(
        flanked_blacklist_temp,
        output_paths.flanked_blacklist.clone(),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("outputs_tests.rs");
}
