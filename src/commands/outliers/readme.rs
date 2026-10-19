//! Human-readable explanations shipped alongside each run's outlier tables.

use super::outputs::{OutputMetadata, OutputPaths};
use crate::shared::io::FinalOutputFiles;
use anyhow::{Context, Result};
use std::{
    fs::File,
    io::{BufWriter, Write},
};

/// Write the file guide and the settings needed to interpret this run's preservation decisions.
///
/// Explanations belong here rather than in long setting-like comments repeated in every table.
/// Register this file with the other outputs so failed runs do not publish a partial guide.
pub(super) fn write_output_readme(
    metadata: &OutputMetadata<'_>,
    paths: &OutputPaths,
    final_outputs: &mut FinalOutputFiles,
) -> Result<()> {
    let temporary = final_outputs.temp_path_for(&paths.readme)?;
    let mut writer = BufWriter::new(File::create(&temporary).context("creating outlier README")?);
    let config = metadata.config;
    writeln!(writer, "cfDNAlab outlier results\n")?;
    writeln!(writer, "Model: {}", config.model.as_str())?;
    writeln!(writer, "Correction target: {}", config.target.as_str())?;
    writeln!(
        writer,
        "Analysed chromosomes: {}",
        metadata.chromosomes.join(", ")
    )?;
    writeln!(
        writer,
        "Eligible positions: {}",
        metadata.eligible_positions
    )?;
    writeln!(
        writer,
        "Tail probability: {:.6e}",
        metadata.tail_probability
    )?;
    writeln!(
        writer,
        "Preservation window: {} bp",
        config.protection_window_size
    )?;
    writeln!(
        writer,
        "Minimum window denominator: {} bases",
        config.min_protection_positions
    )?;
    writeln!(
        writer,
        "Preservation cutoff: strictly greater than {}%",
        config.max_outlier_percentage
    )?;
    writeln!(
        writer,
        "Action above cutoff: {}",
        config.ambiguous_region_action.as_str()
    )?;
    writeln!(
        writer,
        "Flank on each side: {} bp\n",
        metadata.blacklist_flank
    )?;
    writeln!(
        writer,
        "Preserved-region plots requested: {}",
        config.plot_preserved_regions
    )?;
    writeln!(
        writer,
        "Plot context: {} times the grouped span, at least {} bp\n",
        config.preserved_plot_context_multiplier, config.preserved_plot_min_span
    )?;
    writer.write_all(OUTPUT_GUIDE.as_bytes())?;
    writer.flush()?;
    drop(writer);
    final_outputs.record(temporary, paths.readme.clone())?;
    Ok(())
}

const OUTPUT_GUIDE: &str = r#"Where to start
--------------
Read outliers.summary.tsv for the amount of sequence called. Inspect outliers.audit.tsv
for individual calls and why they were selected for correction or left unchanged.
Use outliers.keep_weights.tsv with fcoverage --outlier-weights to apply the correction.
Compare raw and weighted fcoverage tracks to assess the effect on coverage.

The detector targets smaller sample-specific pileups remaining after ordinary blacklisting.
Preservation is a heuristic to avoid correcting broad or clustered elevations. It does not
identify their biological or technical cause. Small biological elevations may be corrected,
and dense technical pileups may be preserved.

Files and usage
---------------
An optional output prefix is followed by a dot before each filename below.

outliers.keep_weights.tsv
    Correction input for fcoverage --outlier-weights. Each row is an exact, unflanked
    call selected for correction. keep_weight is the fraction of a fragment's contribution
    to retain, between 0 and 1. Unlisted positions have implicit weight 1.
    Any overlap applies the full regional weight to the whole fragment, including its
    contribution outside the call. Multiple overlaps use the smallest weight.
    A zero weight excludes the fragment before GC correction, normalization, and counting.
    Use the same fragment filters and blacklists for detection and fcoverage.
    GC correction weights multiply with the outlier weights when both are supplied.

outliers.exact.bed
    Exact, unflanked calls selected for correction, as three-column BED intervals.
    This is an exclusion alternative to fractional weights. A receiving command's
    blacklist rules determine what is excluded. Preserved calls are absent.

outliers.flanked.bed
    Calls extended on both sides, clipped to chromosome bounds, then merged where
    extensions overlap or touch. Preserved merged regions are absent.
    This blacklist also excludes nearby sequence. Only this BED adds flanking to
    the output exclusion coordinates. Weights and the exact BED remain unflanked.

outliers.audit.tsv
    Every exact candidate, including preserved calls. See the column guide below.
    This is an inspection file, not correction input. proposed_keep_weight is the
    proposed correction even when action=preserve means it will not be applied.

outliers.summary.tsv
    Counts and genomic extents for each chromosome and for all analysed chromosomes.
    See the column guide below. These values do not measure removed coverage.

outliers.preserved_regions.tsv
    Always written. Preserved calls grouped into approximate regions for inspection,
    with call counts, exact called bases, selected window bounds, and plot coordinates.
    This table has a header and no data rows when no calls are preserved.
    See the preserved-region column guide below.

outliers.preserved_regions.region_000001.png (and one numbered PNG per group)
    Written only with --plot-preserved-regions. Show surrounding raw coverage and a
    closer view, with the fitted background mean, outlier cutoff, preserved calls,
    and blacklisted gaps. plot_file in the table identifies each figure.
    Coverage uses the same filters as detection. Indexed BAM reads are restricted
    to the recorded plot contexts and do not contribute to detection statistics.
    Each reference bin shows mean and peak coverage among its unblacklisted bases.
    Genuine zeros are included. Grey masked spans are unavailable, not zero coverage.
    The displayed bin width is chosen to keep at most 2,000 bins in the context.

outliers.models.tsv
    Coverage model and calling cutoff for each consecutive --stride interval.
    start/end describe where the model is applied. context_start/context_end describe
    the broader region used to estimate it. A call spanning several intervals can use
    several cutoffs and target means. See the model column guide below.

outliers.histograms.tsv.zst
    Zstandard-compressed raw coverage counts, before weighting or GC correction.
    Each start/end interval is a --stride interval. coverage is an integer fragment
    coverage. observed_positions counts unblacklisted bases at that coverage,
    including zero coverage. Missing coverage values have zero positions.
    fragments_by_midpoint counts accepted fragments whose midpoint is in this interval
    and is not blacklisted. Each fragment is assigned once by midpoint. The count
    repeats for each coverage row within the interval. Do not sum these repetitions.
    A fully blacklisted interval has a single row with zero observed positions.

outliers.global_fit.tsv
    The complete sample-wide coverage histogram and both fitted distributions.
    observed_positions counts bases at each coverage. initial_fitted_positions and
    underlying_fitted_positions give expected counts from the first and second fits.
    Parameter and diagnostic values in the leading comments describe the global fit.
    Local cutoffs in outliers.models.tsv determine calls, except where a local row
    explicitly uses the global fallback.

outliers.global_fit.png
    Global observed and fitted distributions, plus chromosome-position tracks of
    local means and cutoffs. The complete observed tail is retained.

outliers.README.txt
    This guide, with the preservation settings actually used in this run.

Coordinates and comments
------------------------
All intervals are 0-based and half-open, as in BED. Their length is end - start.
For example, [100, 110) contains 10 bases. Eligible positions are unblacklisted bases
on the analysed chromosomes, including bases with zero coverage.
Text tables begin with # key=value comments recording settings or fitted values.
Ignore lines beginning with # when reading the data header and rows.
analysed_chromosomes records the analysed names and reference lengths as JSON objects
with name and size fields. It also identifies analysed chromosomes with no calls.
blacklists records the supplied blacklist paths for reference.

Summary columns
---------------
scope and chromosome distinguish chromosome rows from the total row.
eligible_positions is the number of unblacklisted bases represented by that row.
candidate_calls counts all exact connected intervals that crossed the final cutoff.
actionable_calls counts those selected for correction. preserved_calls counts those
left unchanged. Candidate counts equal actionable plus preserved counts.

candidate_bases, actionable_bases, and preserved_bases sum exact call lengths.
These counts exclude flanks and the gaps between calls. The corresponding percentage
columns divide by eligible_positions in the same row. A zero denominator gives NaN.

longest_candidate_call_bp, longest_actionable_call_bp, and longest_preserved_call_bp
give the longest SINGLE exact call in each category, before flanking. They are lengths
in bases, not call counts, total called bases, window sizes, or merged-region lengths.
The total row takes the maximum over chromosomes. It does not sum these lengths.

actionable_merged_regions and preserved_merged_regions count the regions formed after
flanking and merging all candidates. Many exact calls can belong to one merged region.
actionable_flanked_span sums full merged spans selected for correction, including
blacklisted bases and gaps between calls. actionable_flanked_eligible_bases counts
only their unblacklisted bases. actionable_flanked_eligible_percentage divides that
count by eligible_positions in the row.

model_cores counts --stride intervals. global_fallback_cores counts those using the
sample-wide model. poisson_boundary_cores counts those whose second-fit dispersion
is zero, including all ZIP fits. A ZINB fit can also select this boundary.

Preserved-region columns
------------------------
region_id identifies the reporting group. chromosome/start/end bound its exact calls,
including the gaps between them but without adding flanks. These are approximate
locations of preserved signal, not inferred biological breakpoints.
calls counts member exact calls. called_bases sums only their exact lengths.
maximum_window_percentage is the largest selected percentage among member calls.
It counts all candidate bases in that window, including candidates outside this group.
window_start/window_end bound their selected preservation windows.

Each existing merged preserved region stays together. Regions whose selected window
spans overlap are joined transitively for reporting. Merely touching window spans do
not join groups. Broader plot contexts never join groups or change corrections.

plot_start/plot_end record the wider context, even without plots. It is centered on
the grouped span, widened by --preserved-plot-context-multiplier, and has at least
--preserved-plot-min-span bases where chromosome bounds allow. It also includes the
selected preservation windows. plot_file is the PNG filename, or NA when not requested.
The close-up panel adds a small margin around the exact-call span, within that context.

Audit columns and preservation
------------------------------
start/end give the exact candidate, before flanking.
observed_mass sums raw fragment coverage over that exact call. target_mass sums
the requested target coverage over the same bases. Their units are fragment-bases.
mean_coverage is observed_mass divided by call length. max_coverage is the highest
coverage within this exact call, not the surrounding model context or processing tile.
proposed_keep_weight is target_mass / observed_mass, capped at 1.

window_start/window_end give the selected physical preservation window.
window_eligible_positions counts its unblacklisted bases, including sampling zeros.
window_called_positions counts exact candidate bases within it, before preservation.
window_score = 100 * window_called_positions /
               max(window_eligible_positions, minimum window denominator).
window_direction is forward or backward. call_flagged is true when the score is
strictly greater than the cutoff printed at the top of this guide.

A forward window starts at every exact call's start. Each call receives the highest
score from overlapping forward windows, including earlier calls' windows. If any of
these windows meets the minimum eligible-base count, ignore all windows below that
minimum when taking the maximum. Thus a short call can be preserved because many
nearby calls collectively exceed the cutoff in a sufficiently supported window.
If no overlapping forward window meets the minimum, also evaluate a window ending
at the call's exact end. Use this backward window if it has enough eligible bases.
Otherwise retain the maximum among the undersupported forward windows with the
denominator floor. Windows stop at chromosome bounds and are never lengthened to
collect bases across blacklisted gaps.

All candidates are then flanked. Overlapping or touching flanks are merged, including
connections through chains of calls. merged_start/merged_end identify the shared
merged region. merged_max_score is the largest selected score among its member calls.
If any member is flagged, the configured action applies to all calls in this region.
Preservation does not trigger another model fit or restore excluded values to the fit.

action is correct or preserve. decision_reason explains it:
    density_above_cutoff: this call's selected window exceeded the cutoff.
    merged_with_preserved_call: this call's window did not exceed the cutoff,
        but its flanks connected it to a call requesting preservation.
    density_at_or_below_cutoff: no member requested preservation.
    flagged_region_correction_requested: a member exceeded the cutoff, but the
        configured action explicitly requested correction anyway.
first_model_core_start/last_model_core_start identify the first and last model rows
overlapped by this call. Intermediate model intervals also apply if present.

Model columns and diagnostics
-----------------------------
model_source is local or global_fallback. The latter means the complete chromosome
could not provide the requested context size and fragment count.
eligible_positions and fragments_by_midpoint describe the recorded context, even for
global fallback. fragments_by_midpoint counts accepted fragments with an unblacklisted
midpoint between context_start and context_end. It does not count every fragment
overlapping the context. Contexts overlap, so do not sum these counts across rows.

initial_* values come from the first fit to the complete histogram.
initial_threshold (T1) excludes coverage >= T1 from the second fit.
second_fit_retained_positions counts positions below T1 used for that second fit.
The second fit accounts for the excluded upper tail. Its underlying_* parameters
describe the full fitted distribution, rather than only values below T1.
final_threshold (T2) is the inclusive calling cutoff: coverage >= T2 is called.

initial_component_mean and underlying_component_mean describe the coverage component.
The model's extra zeros are separate. underlying_mean includes those extra zeros:
    underlying_mean = (1 - underlying_zero_inflation) * underlying_component_mean.
The local-mean correction target uses underlying_component_mean.
The *_zero_inflation values are probabilities assigned to the extra-zero component.
The *_dispersion values allow variance above a Poisson distribution in the coverage
component: variance = mean + dispersion * mean^2. Zero is the Poisson boundary,
recorded in the corresponding *_poisson_boundary flags.

observed_positive_coverage_mean and observed_positive_coverage_variance describe
the original positions with coverage > 0, including extreme values. The variance is
the population mean squared deviation, not uncertainty in an estimated parameter.
underlying_expected_positive_coverage_variance is the variance among positive values
expected from the second-fit model. positive_coverage_variance_ratio divides the
observed variance by that expectation. Near 1 means similar spread. Above 1 means
more observed variation, and below 1 means less. A high ratio can reflect pileups
or a poor description of ordinary coverage. This diagnostic does not change calls.
Observed positive moments and their ratio are NaN if an incomplete context has no
positive coverage and therefore uses the global fallback model.

observed_initial_tail_positions and expected_initial_tail_positions compare the
original count of bases at coverage >= T1 with the first-fit model expectation.
observed_final_tail_positions and expected_final_tail_positions do the same at T2
using the second-fit model. Excess counts can include the pileups being detected.
Neighbouring bases share fragments, so these comparisons are not independent tests.

Genomic extent and coverage reduction
------------------------------------
Called-base percentages measure genomic extent, not the fraction of coverage removed.
Full-fragment weighting also changes coverage outside exact calls. The masses in the
audit do not measure the complete change in a weighted coverage track.
Compare raw and weighted fcoverage outputs with identical filters and blacklists,
initially without GC correction, genomic scaling, or length normalization. Inspect
the largest pileups and the profiles or features that matter for the analysis.
"#;
