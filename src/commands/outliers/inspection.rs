//! Actionable call files and complete audit/summary outputs for BAM inspection.

use super::{
    outliers::CalledRegion,
    outputs::{OutputMetadata, OutputPaths, write_common_metadata},
    preservation::{ChromosomeDecisions, IntervalLengths, evaluate_regions},
    preserved_regions::{PreservedRegion, collect_preserved_regions},
    striding::{CoreHistogram, CoreModel, ModelSource},
};
use crate::shared::{bam::Contigs, interval::Interval, io::FinalOutputFiles};
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

/// Counts use exact eligible called bases. Flanked spans and eligible flanked bases are separate.
#[derive(Debug, Default)]
pub(super) struct CallSummary {
    /// All eligible positions in this summary's chromosome or combined selected chromosomes.
    eligible_positions: u64,
    pub(super) candidate_calls: u64,
    pub(super) actionable_calls: u64,
    pub(super) preserved_calls: u64,
    candidate_bases: u64,
    pub(super) actionable_bases: u64,
    pub(super) preserved_bases: u64,
    longest_candidate: u32,
    longest_actionable: u32,
    longest_preserved: u32,
    actionable_merged_regions: u64,
    preserved_merged_regions: u64,
    /// Physical length of actionable merged flanks, including masked bases and gaps between calls.
    actionable_flanked_span: u64,
    /// Unblacklisted part of the same merged flanks, used for eligible-genome percentages.
    actionable_flanked_eligible_bases: u64,
    model_cores: u64,
    global_fallback_cores: u64,
    poisson_boundary_cores: u64,
}

impl CallSummary {
    fn add(&mut self, other: &Self) {
        self.eligible_positions += other.eligible_positions;
        self.candidate_calls += other.candidate_calls;
        self.actionable_calls += other.actionable_calls;
        self.preserved_calls += other.preserved_calls;
        self.candidate_bases += other.candidate_bases;
        self.actionable_bases += other.actionable_bases;
        self.preserved_bases += other.preserved_bases;
        self.longest_candidate = self.longest_candidate.max(other.longest_candidate);
        self.longest_actionable = self.longest_actionable.max(other.longest_actionable);
        self.longest_preserved = self.longest_preserved.max(other.longest_preserved);
        self.actionable_merged_regions += other.actionable_merged_regions;
        self.preserved_merged_regions += other.preserved_merged_regions;
        self.actionable_flanked_span += other.actionable_flanked_span;
        self.actionable_flanked_eligible_bases += other.actionable_flanked_eligible_bases;
        self.model_cores += other.model_cores;
        self.global_fallback_cores += other.global_fallback_cores;
        self.poisson_boundary_cores += other.poisson_boundary_cores;
    }

    /// Express genomic base counts relative to this summary's full eligible genomic span.
    ///
    /// This denominator differs from a call's preservation-window denominator and does not measure
    /// the coverage reduction caused by weighting complete overlapping fragments.
    pub(super) fn percentage(&self, bases: u64) -> f64 {
        if self.eligible_positions == 0 {
            f64::NAN
        } else {
            100.0 * bases as f64 / self.eligible_positions as f64
        }
    }

    fn write(&self, writer: &mut impl Write, scope: &str, chromosome: &str) -> Result<()> {
        writeln!(
            writer,
            "{}",
            [
                scope.to_string(),
                chromosome.to_string(),
                self.eligible_positions.to_string(),
                self.candidate_calls.to_string(),
                self.actionable_calls.to_string(),
                self.preserved_calls.to_string(),
                self.candidate_bases.to_string(),
                self.percentage(self.candidate_bases).to_string(),
                self.actionable_bases.to_string(),
                self.percentage(self.actionable_bases).to_string(),
                self.preserved_bases.to_string(),
                self.percentage(self.preserved_bases).to_string(),
                self.longest_candidate.to_string(),
                self.longest_actionable.to_string(),
                self.longest_preserved.to_string(),
                self.actionable_merged_regions.to_string(),
                self.preserved_merged_regions.to_string(),
                self.actionable_flanked_span.to_string(),
                self.actionable_flanked_eligible_bases.to_string(),
                self.percentage(self.actionable_flanked_eligible_bases)
                    .to_string(),
                self.model_cores.to_string(),
                self.global_fallback_cores.to_string(),
                self.poisson_boundary_cores.to_string()
            ]
            .join("\t")
        )?;
        Ok(())
    }
}

/// Summarize each exact call once and each merged flank once, independently of fragment counters.
fn summarize_calls(
    regions: &[CalledRegion],
    decisions: &ChromosomeDecisions,
    eligible_positions: u64,
    masked: &IntervalLengths,
    models: &[CoreModel],
) -> CallSummary {
    let mut summary = CallSummary {
        eligible_positions,
        model_cores: models.len() as u64,
        global_fallback_cores: models
            .iter()
            .filter(|model| model.source == ModelSource::GlobalFallback)
            .count() as u64,
        poisson_boundary_cores: models
            .iter()
            .filter(|model| model.model.underlying_fit.dispersion == 0.0)
            .count() as u64,
        ..CallSummary::default()
    };
    for (region, decision) in regions.iter().zip(&decisions.calls) {
        let length = region.interval.len();
        summary.candidate_calls += 1;
        summary.candidate_bases += length as u64;
        summary.longest_candidate = summary.longest_candidate.max(length);
        if decisions.merged_regions[decision.merged_region_index].preserved {
            summary.preserved_calls += 1;
            summary.preserved_bases += length as u64;
            summary.longest_preserved = summary.longest_preserved.max(length);
        } else {
            summary.actionable_calls += 1;
            summary.actionable_bases += length as u64;
            summary.longest_actionable = summary.longest_actionable.max(length);
        }
    }
    for merged in &decisions.merged_regions {
        if merged.preserved {
            summary.preserved_merged_regions += 1;
        } else {
            summary.actionable_merged_regions += 1;
            summary.actionable_flanked_span += merged.interval.len() as u64;
            summary.actionable_flanked_eligible_bases +=
                merged.interval.len() as u64 - masked.within(merged.interval);
        }
    }
    summary
}

/// Write every candidate's evidence and action. Model-core starts identify all relevant fit rows.
fn write_audit_rows(
    writer: &mut impl Write,
    chromosome: &str,
    regions: &[CalledRegion],
    decisions: &ChromosomeDecisions,
    metadata: &OutputMetadata<'_>,
) -> Result<()> {
    for (region, decision) in regions.iter().zip(&decisions.calls) {
        let merged = &decisions.merged_regions[decision.merged_region_index];
        let reason = match (merged.preserved, merged.flagged, decision.flagged) {
            (true, _, true) => "density_above_cutoff",
            (true, _, false) => "merged_with_preserved_call",
            (false, true, _) => "flagged_region_correction_requested",
            _ => "density_at_or_below_cutoff",
        };
        let window = decision.window;
        // Mean/maximum coverage and both masses describe this exact call. The selected window
        // and merged flank span are separate audit fields used to explain its preservation action
        writeln!(
            writer,
            "{}",
            [
                chromosome.to_string(),
                region.interval.start().to_string(),
                region.interval.end().to_string(),
                region.observed_mass.to_string(),
                region.target_mass.to_string(),
                (region.observed_mass / region.interval.len() as f64).to_string(),
                region.maximum_coverage.to_string(),
                region.keep_weight().to_string(),
                window.score.to_string(),
                window.interval.start().to_string(),
                window.interval.end().to_string(),
                window.eligible_positions.to_string(),
                window.called_positions.to_string(),
                if window.backward {
                    "backward"
                } else {
                    "forward"
                }
                .to_string(),
                decision.flagged.to_string(),
                merged.interval.start().to_string(),
                merged.interval.end().to_string(),
                merged.maximum_score.to_string(),
                if merged.preserved {
                    "preserve"
                } else {
                    "correct"
                }
                .to_string(),
                reason.to_string(),
                (region.interval.start() / metadata.config.stride * metadata.config.stride)
                    .to_string(),
                ((region.interval.end() - 1) / metadata.config.stride * metadata.config.stride)
                    .to_string()
            ]
            .join("\t")
        )?;
    }
    Ok(())
}

/// Publish actionable exact weights/BEDs together with all-candidate evidence and genomic totals.
pub(super) fn write_call_outputs(
    metadata: &OutputMetadata<'_>,
    contigs: &Contigs,
    regions_by_chromosome: &FxHashMap<String, Vec<CalledRegion>>,
    blacklists: &FxHashMap<String, Vec<Interval<u64>>>,
    histograms: &FxHashMap<String, Vec<CoreHistogram>>,
    models: &FxHashMap<String, Vec<CoreModel>>,
    paths: &OutputPaths,
    final_outputs: &mut FinalOutputFiles,
) -> Result<(CallSummary, Vec<PreservedRegion>)> {
    let destinations = [
        &paths.keep_weights,
        &paths.exact_blacklist,
        &paths.flanked_blacklist,
        &paths.audit,
        &paths.summary,
        &paths.preserved_regions,
    ];
    let temporary = destinations
        .iter()
        .map(|path| final_outputs.temp_path_for(path))
        .collect::<Result<Vec<_>>>()?;
    let mut writers = temporary
        .iter()
        .map(|path| Ok(BufWriter::new(File::create(path)?)))
        .collect::<Result<Vec<_>>>()?;
    let [
        keep,
        exact,
        flanked,
        audit,
        summary_writer,
        preserved_writer,
    ] = writers.as_mut_slice()
    else {
        unreachable!()
    };
    for writer in [
        &mut *keep,
        &mut *exact,
        &mut *flanked,
        &mut *audit,
        &mut *summary_writer,
        &mut *preserved_writer,
    ] {
        write_common_metadata(writer, metadata)?;
    }
    writeln!(
        keep,
        "# omitted_keep_weight=1.0\n# fragment_overlap_rule=minimum_keep_weight\nchromosome\tstart\tend\tkeep_weight"
    )?;
    writeln!(exact, "# output=actionable_exact_called_regions")?;
    writeln!(flanked, "# output=actionable_merged_flanked_regions")?;
    writeln!(
        audit,
        "chromosome\tstart\tend\tobserved_mass\ttarget_mass\tmean_coverage\tmax_coverage\tproposed_keep_weight\twindow_score\twindow_start\twindow_end\twindow_eligible_positions\twindow_called_positions\twindow_direction\tcall_flagged\tmerged_start\tmerged_end\tmerged_max_score\taction\tdecision_reason\tfirst_model_core_start\tlast_model_core_start"
    )?;
    writeln!(
        summary_writer,
        "scope\tchromosome\teligible_positions\tcandidate_calls\tactionable_calls\tpreserved_calls\tcandidate_bases\tcandidate_percentage\tactionable_bases\tactionable_percentage\tpreserved_bases\tpreserved_percentage\tlongest_candidate_call_bp\tlongest_actionable_call_bp\tlongest_preserved_call_bp\tactionable_merged_regions\tpreserved_merged_regions\tactionable_flanked_span\tactionable_flanked_eligible_bases\tactionable_flanked_eligible_percentage\tmodel_cores\tglobal_fallback_cores\tpoisson_boundary_cores"
    )?;
    writeln!(
        preserved_writer,
        "region_id\tchromosome\tstart\tend\tcalls\tcalled_bases\tmaximum_window_percentage\twindow_start\twindow_end\tplot_start\tplot_end\tplot_file"
    )?;
    let mut total = CallSummary::default();
    let mut preserved_regions = Vec::new();
    for chromosome in metadata.chromosomes {
        let regions = regions_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing calls for {chromosome}"))?;
        let chromosome_length = contigs
            .contigs
            .get(chromosome)
            .with_context(|| format!("missing contig {chromosome}"))?
            .1;
        let blacklist = blacklists
            .get(chromosome)
            .map(Vec::as_slice)
            .unwrap_or_default();
        let decisions = evaluate_regions(metadata.config, regions, blacklist, chromosome_length)
            .with_context(|| format!("preservation decisions for {chromosome}"))?;
        write_audit_rows(audit, chromosome, regions, &decisions, metadata)?;
        let reporting_groups = collect_preserved_regions(
            chromosome,
            regions,
            &decisions,
            metadata.config,
            chromosome_length,
            preserved_regions.len() + 1,
        )?;
        for group in &reporting_groups {
            let plot_file = if metadata.config.plot_preserved_regions {
                paths
                    .preserved_plot(&group.region_id)
                    .file_name()
                    .context("preserved plot path has no filename")?
                    .to_string_lossy()
                    .into_owned()
            } else {
                "NA".to_string()
            };
            writeln!(
                preserved_writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                group.region_id,
                chromosome,
                group.interval.start(),
                group.interval.end(),
                group.call_indices.len(),
                group.called_bases,
                group.maximum_window_percentage,
                group.window_span.start(),
                group.window_span.end(),
                group.plot_interval.start(),
                group.plot_interval.end(),
                plot_file
            )?;
        }
        preserved_regions.extend(reporting_groups);
        for (region, decision) in regions.iter().zip(&decisions.calls) {
            if decisions.merged_regions[decision.merged_region_index].preserved {
                continue;
            }
            writeln!(
                keep,
                "{}\t{}\t{}\t{}",
                chromosome,
                region.interval.start(),
                region.interval.end(),
                region.keep_weight()
            )?;
            writeln!(
                exact,
                "{}\t{}\t{}",
                chromosome,
                region.interval.start(),
                region.interval.end()
            )?;
        }
        for merged in &decisions.merged_regions {
            if !merged.preserved {
                writeln!(
                    flanked,
                    "{}\t{}\t{}",
                    chromosome,
                    merged.interval.start(),
                    merged.interval.end()
                )?;
            }
        }
        let chromosome_histograms = histograms
            .get(chromosome)
            .with_context(|| format!("missing histograms for {chromosome}"))?;
        // Model cores partition the chromosome, so their eligible-position counts include each
        // eligible position once, even though the broader fitting contexts overlap
        let eligible_positions =
            chromosome_histograms
                .iter()
                .try_fold(0_u64, |total, histogram| -> Result<u64> {
                    total
                        .checked_add(histogram.eligible_positions()?)
                        .context("chromosome eligible-position overflow")
                })?;
        let chromosome_models = models
            .get(chromosome)
            .with_context(|| format!("missing models for {chromosome}"))?;
        let summary = summarize_calls(
            regions,
            &decisions,
            eligible_positions,
            &IntervalLengths::new(blacklist.to_vec()),
            chromosome_models,
        );
        summary.write(summary_writer, "chromosome", chromosome)?;
        total.add(&summary);
    }
    total.write(summary_writer, "total", "all")?;
    for writer in &mut writers {
        writer.flush()?;
    }
    drop(writers);
    for (temporary, destination) in temporary.into_iter().zip(destinations) {
        final_outputs.record(temporary, destination.clone())?;
    }
    Ok((total, preserved_regions))
}

#[cfg(test)]
mod tests {
    include!("inspection_tests.rs");
}
