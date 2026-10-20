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
        "Minimum eligible positions for window percentages: {}",
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

const OUTPUT_GUIDE: &str = r#"Start here
----------
This run looked for unusually high fragment coverage.
Coverage is the number of fragments covering a base. A call is a stretch of adjacent
bases whose coverage reached the outlier cutoff. A flank is an extra buffer added
on each side of a call in the flanked blacklist.

Calls are divided into two groups:
- Actionable calls are selected for correction.
- Preserved calls are left unchanged because they may belong to a broader elevation.
  This does not establish whether the elevation is biological or technical.

Read `outliers.summary.tsv` to see how much sequence was called. Inspect individual
calls in `outliers.audit.tsv`. Use `outliers.keep_weights.tsv` with
`fcoverage --outlier-weights` to apply the correction. Compare coverage before and
after weighting to see its effect on your analysis.

The command targets extremely high coverage over short genomic intervals after ordinary
blacklisting. It cannot identify the cause from coverage alone. Small biological
events may be corrected, and broad technical pileups may be preserved.

Files and usage
---------------
If you supplied an output prefix, it appears before each filename below.

`outliers.keep_weights.tsv`
    Input for `fcoverage --outlier-weights`. Each row gives an exact call selected
    for correction and its `keep_weight`, the fraction of a fragment's contribution
    to retain. By default, it aims to reduce regional coverage toward the estimated
    background mean. For example, 0.1 retains 10%. A zero weight excludes the fragment
    before GC correction, normalization calculations, and counting.
    Even a 1 bp overlap applies the full weight to the whole fragment, including
    coverage outside the call. Multiple overlaps use the smallest weight.
    Fragments without an overlap are unchanged. No flanks are added to these calls.
    Use the same fragment filters and blacklists for detection and `fcoverage`.
    When GC correction is also supplied, its weights multiply the outlier weights.

`outliers.exact.bed`
    Calls selected for correction, without flanks. Use this as a blacklist to
    exclude called regions instead of applying fractional weights. The command
    receiving the blacklist determines whether it excludes bases, fragments,
    or another counting unit.

`outliers.flanked.bed`
    The blacklist alternative with a buffer added on both sides of each call.
    Buffers stop at chromosome ends. Overlapping or touching buffers are joined.
    Preserved regions are omitted from both BED files and the keep-weight table.

`outliers.audit.tsv`
    Every call, including preserved calls, with its coverage, proposed correction,
    and final decision. Use this to investigate why a call was corrected or preserved.
    This file is for inspection. It is not an input for applying corrections.

`outliers.summary.tsv`
    Counts of calls and bases for each chromosome and for the complete run.
    Exact called bases and the extent of the flanked blacklist are reported separately.
    These percentages describe genomic sequence, not how much coverage is removed.

`outliers.preserved_regions.tsv`
    Preserved calls grouped into approximate regions for inspection. Each row gives
    the group's location, number of calls, called bases, and plotting coordinates.
    The table is always written. It has no data rows if no calls were preserved.

`outliers.preserved_regions.region_000001.png` (one numbered PNG per group)
    Added with `--plot-preserved-regions`. Each figure shows surrounding coverage and
    a separately binned closer view, with the active-component mean, outlier cutoff,
    and all nearby candidates colored by preservation or correction.
    Blue segments show mean coverage within bins. Grey marks show the bin peaks so
    narrow pileups remain visible. The bin width is printed above each panel.
    Blacklisted spans are shaded grey. They are excluded from the mean and peak,
    while genuine zero-coverage bases are included.
    Plotting reads the BAM again around these groups using the detection filters.
    These extra reads are not added to the detection statistics. Progress logs identify
    each group and its requested span. Many overlapping contexts can cause substantial rereading.
    Reusing a prefix does not remove old optional PNGs. Use only files listed in this run's
    preserved-region table, or use a fresh prefix for each attempt.

`outliers.models.tsv`
    The coverage model and cutoff used in each consecutive `--stride` interval.
    Each model is usually estimated from a larger surrounding region, its context.
    Use this table to find the background mean and cutoff at a particular location.
    The model-column reference below also explains the fit diagnostics.

`outliers.global_fit.png`
    A quick overview of the observed coverage distribution and fitted models.
    The first two panels compare the sample-wide fit with observed coverage when a
    global fallback is required. Otherwise they show observations and local threshold
    summaries, explicitly labeled "global fit not needed".
    The bottom panel shows the fitted means and applied cutoffs along labeled chromosomes.
    Its x coordinate is the cumulative model-core index in selected chromosome order,
    matching consecutive rows in `outliers.models.tsv`. The unconditional mean is diagnostic.
    The active-component mean is the correction target only for `--target local-mean`.
    Sample-wide cutoffs in the first two panels are references, not universal cutoffs
    applied everywhere. Local cutoffs are recorded in `outliers.models.tsv`.

`outliers.global_fit.tsv`
    Observed counts of bases at each coverage level across the analysed chromosomes,
    alongside counts expected from the first and second sample-wide fits.
    `coverage` is the number of fragments covering a base. `observed_positions` counts
    bases with that coverage. `initial_fitted_positions` and `underlying_fitted_positions`
    give the expected counts from the first and second fits, respectively.
    A global model is fitted only if a chromosome needs fallback. The leading
    `global_fit_status=not_needed` comment identifies skipped fits. Then expected-count columns are `NA`
    and global fit parameters are omitted. Observed counts are always retained.

`outliers.histograms.tsv.zst`
    Zstandard-compressed counts of bases at each raw coverage level, split by
    `--stride` interval. Use these for detailed checks without rereading the BAM.
    `chromosome`, `start`, and `end` locate the interval. `coverage` gives a coverage
    level, and `observed_positions` counts its eligible positions at that level.
    Zero-coverage bases are included. Omitted coverage levels have no observed bases.
    `fragments_by_midpoint` counts accepted fragments whose midpoint falls in this
    interval and is not blacklisted. This count repeats on every row for the interval.
    Do not add those repetitions together. A fully blacklisted interval has a single
    row with zero `observed_positions`.

`outliers.README.txt`
    This file, including the settings used in this run.

Coordinates and table comments
------------------------------
`start` is 0-based and `end` is excluded, as in BED files. The length is `end - start`.
For example, [100, 110) contains 10 bases. Eligible positions are positions outside
the supplied blacklists on the analysed chromosomes. They include positions with zero coverage.

Lines starting with `#` record settings or fitted values. Skip them when reading
a table's header and rows. `analysed_chromosomes` records the analysed chromosome
names and lengths, including chromosomes with no calls. `blacklists` records the
supplied blacklist paths.

Summary columns
---------------
`scope` is `chromosome` for a chromosome row and `total` for the whole run.
`chromosome` names the chromosome, or is `all` for the total row.
`eligible_positions` counts the eligible positions in that row, including those with
zero coverage.

`candidate_calls` counts all calls before deciding which to preserve.
`actionable_calls` counts calls selected for correction. `preserved_calls` counts
calls left unchanged. Together they add up to `candidate_calls`.

`candidate_bases`, `actionable_bases`, and `preserved_bases` count the bases in those
calls. They exclude flanks and gaps between calls. `candidate_percentage`,
`actionable_percentage`, and `preserved_percentage` express those counts as percentages
of `eligible_positions` in the same row. `NaN` means the row has no eligible bases.

`longest_candidate_call_bp`, `longest_actionable_call_bp`, and `longest_preserved_call_bp`
give the longest individual call in each group, in bases and before flanking.
The total row reports the longest call across chromosomes, rather than adding lengths.

`actionable_merged_regions` and `preserved_merged_regions` count groups formed by
adding flanks and joining overlapping or touching flanks. A group can contain many calls.
`actionable_flanked_span` is the total length of groups selected for the flanked
blacklist, including blacklisted bases and gaps between their calls.
`actionable_flanked_eligible_bases` counts only their eligible positions.
`actionable_flanked_eligible_percentage` expresses that count as a percentage of
`eligible_positions` in the same row.

`model_cores` counts the `--stride` intervals that received a coverage model.
`global_fallback_cores` counts intervals that used the sample-wide model because their
chromosome could not provide the requested span and fragment count for a local fit.
`poisson_boundary_cores` counts intervals where the second fit found no extra
negative-binomial variation. All ZIP fits fall in this group. ZINB can also do so.

Preserved-region columns
------------------------
`region_id` identifies a group of preserved calls. `chromosome`, `start`, and `end`
give its location, including gaps between calls but without adding flanks. These bounds
show approximately where preserved signal lies. They do not estimate biological breakpoints.
`calls` counts the group's calls. `called_bases` counts only bases in those calls.

Calls joined by their flanks stay in the same inspection group. Groups are then joined
when the spans containing their selected preservation windows overlap. This can connect
several groups through a chain of overlaps. Windows that merely touch do not join groups.
Grouping for inspection does not change any correction decisions.

`maximum_window_percentage` is the highest selected window score among these calls.
That score counts all candidate bases in its window, including calls outside this group.
`window_start` and `window_end` bound all the selected windows for the group's calls.
They may describe a wider span than the single window giving the highest score.

`plot_start` and `plot_end` give the wider context for a coverage plot, even when plots
were not requested. Its width is controlled by `--preserved-plot-context-multiplier`
and `--preserved-plot-min-span`. It also includes the selected preservation windows
and stops at chromosome ends. `plot_file` gives the PNG filename, or `NA` without plots.

Audit columns
-------------
`chromosome`, `start`, and `end` locate an individual call before flanking.
`observed_mass` is the sum of raw coverage across its bases. For example, 100 bases
with coverage 20 have mass 2,000. This is not a count of distinct fragments.
`target_mass` sums the requested target coverage over the same bases.
`mean_coverage` is `observed_mass` divided by call length. `max_coverage` is the highest
coverage within this call, not within the larger fitting context or processing tile.
`proposed_keep_weight` is `target_mass / observed_mass`, capped at 1. It is recorded
even for preserved calls, where the proposed correction is not applied.

`window_start` and `window_end` locate the window selected for this call's decision.
`window_eligible_positions` counts eligible positions in that window, including those
with zero coverage.
`window_called_positions` counts candidate bases inside it, before preservation.
`window_score` is the called-base percentage used for the decision:
    100 * `window_called_positions` /
          max(`window_eligible_positions`, minimum eligible positions printed above)
`window_direction` is `forward` or `backward`. `call_flagged` is `true` when
`window_score` is strictly above the preservation cutoff printed above.
A short call can be preserved because nearby calls collectively exceed this cutoff,
or because its flanks connect it to a preserved call.

A forward window starts at each call's start. Each call considers its own window and
all other forward windows that overlap it. If any has enough eligible positions,
use the highest percentage from those windows. Ignore windows below the minimum.
If none has enough, also check a backward window ending at the call's end. Use it
if it has enough eligible positions. Otherwise keep the highest forward-window score,
using the minimum denominator in the formula above. Windows stop at chromosome ends.
Blacklisted bases are skipped without extending the window to collect more bases.

All calls then receive flanks. Overlapping or touching flanks are joined.
`merged_start` and `merged_end` bound the joined region containing this call.
`merged_max_score` is the highest selected score among calls in that region.
If any call is flagged, the configured action applies to all calls in the joined region.
Preservation does not change the background fit or its cutoff.

`action` is `correct` or `preserve`. `decision_reason` explains the action:
- `density_above_cutoff`: this call's selected window exceeded the cutoff.
- `merged_with_preserved_call`: this call's window did not exceed the cutoff, but its
  flanks connected it to a preserved call.
- `density_at_or_below_cutoff`: no call in the joined region exceeded the cutoff.
- `flagged_region_correction_requested`: a call exceeded the cutoff, but
  `--ambiguous-region-action correct` requested correction anyway.

`first_model_core_start` and `last_model_core_start` identify the first and last
`start` values in `outliers.models.tsv` used within this call. Any model intervals
between them also apply. A long call can therefore have more than one background mean
or cutoff contributing to its proposed weight.

Model columns
-------------
ZINB means zero-inflated negative binomial. It allows more coverage variation than a
Poisson model. ZIP means zero-inflated Poisson and is the simpler comparison model.
Both account for additional zero-coverage positions, as explained below.

`chromosome`, `start`, and `end` locate the interval where this model was applied.
`context_start` and `context_end` locate the larger region used to estimate it.
`model_source` is `local` when that chromosome supplied enough context. It is
`global_fallback` when the whole chromosome could not meet the requested context span
and fragment count, so the sample-wide model was used.

`eligible_positions` and `fragments_by_midpoint` describe the recorded chromosome
context, including for fallback rows. `fragments_by_midpoint` counts accepted fragments
whose midpoint is between `context_start` and `context_end` and outside the supplied
blacklists. It does not count all fragments overlapping the context. Contexts overlap,
so do not sum these counts across model rows.

The model is fitted twice. The first fit estimates an initial cutoff, T1. Coverage
at or above T1 is omitted from the second fit so extreme values have less influence
on the estimated background. The second fit accounts for that omission when estimating
the full coverage distribution. It supplies the final calling cutoff, T2.

`initial_threshold` is T1. `second_fit_retained_positions` counts bases below T1 used
for the second fit. `final_threshold` is T2. A base is called when its raw coverage
is at or above T2. Columns starting with `initial_` describe the first fit, while
those starting with `underlying_` describe the background estimated by the second fit.

Both models separate extra zeros from a coverage component that can itself produce
zeros by ordinary sampling. Extra zeros are a fitted description of the distribution.
They do not identify specific unmappable bases or reference gaps.
`initial_component_mean` and `underlying_component_mean` are the coverage-component
means from the two fits. The `local-mean` correction target uses
`underlying_component_mean`. `underlying_mean` also includes the extra zeros:
    `underlying_mean` = (1 - `underlying_zero_inflation`) * `underlying_component_mean`
`initial_zero_inflation` and `underlying_zero_inflation` give the fitted probabilities
of belonging to the extra-zero component.

Model fit diagnostics
---------------------
These columns help assess whether the model describes ordinary coverage well enough.
They do not change the calls or prove that a called region is a technical artifact.

`initial_dispersion` and `underlying_dispersion` describe the extra variation allowed
in the coverage component of each fit. Its variance is mean + dispersion * mean^2.
A dispersion of zero gives the Poisson model. `initial_poisson_boundary` and
`underlying_poisson_boundary` record whether that happened. A ZINB fit may select
this simpler model when extra variation is not supported.

`observed_positive_coverage_mean` and `observed_positive_coverage_variance` summarize
bases with coverage above zero, including extreme pileups. Variance measures how much
the coverage values differ, not uncertainty in the fitted mean.
`underlying_expected_positive_coverage_variance` gives the variance expected among
positive values under the second-fit model. `positive_coverage_variance_ratio` divides
the observed variance by that expectation. Near 1 means similar spread, above 1 means
more observed variation, and below 1 means less. A high ratio can reflect either
pileups or a poor fit to ordinary coverage. The observed mean, observed variance, and
ratio are `NaN` when an all-zero context requires the sample-wide fallback.

`observed_initial_tail_positions` counts bases at or above T1.
`expected_initial_tail_positions` gives the count expected from the first fit.
`observed_final_tail_positions` and `expected_final_tail_positions` make the same
comparison at T2 using the second fit. More calls than expected can reflect pileups.
A widespread excess across contexts can also suggest that the model's cutoff is too
low for ordinary coverage. Neighbouring bases share fragments, so these comparisons
are not independent statistical tests.

`final_threshold_exceeds_initial` flags T2 > T1. Correcting for truncation can raise
rather than lower the fitted threshold. This is not itself an error and is never clamped.
`initial_optimizer_diagnostics` and `underlying_optimizer_diagnostics` are JSON objects
for ZINB, or `NA` for ZIP. They record converged and unfinished starts and ranges across
valid converged candidates within the recorded average negative-log-likelihood tolerance
of the selected fit, including the explicit ZIP boundary when comparable. Ranges cover
component mean, dispersion, structural-zero probability, and that stage's tail threshold.
They describe numerical agreement, not confidence intervals or proof of identification.

Runtime, filters, and reference gaps
-----------------------------------
Detection scans the selected BAM regions twice. Histograms and model fitting follow the
first pass. Calling and preservation follow the second. Optional preserved-region plots
reread their local contexts, in bounded chunks, with the same fragment filters.
Unmapped, secondary, supplementary, duplicate, and QC-failed alignments are excluded.
MAPQ and fragment length limits apply. Paired detection also requires mapped mates on
the same chromosome with opposite strands. Proper-pair flags are required only when
requested. With `--reads-are-fragments`, each accepted alignment is counted separately.

The BAM supplies reference lengths, but the detector does not read reference sequence.
Unmasked reference gaps remain eligible zero-coverage positions. Alignment D/N gaps
contribute no coverage. `--ignore-gap` also excludes the inter-mate gap from coverage.
These gaps remain eligible unless blacklisted. Blacklisting masks positions, not whole
fragments: fragments crossing a mask can still contribute outside it. `fcoverage` likewise
masks positional coverage, whereas outlier weights apply to entire overlapping fragments.

A short or unsupported chromosome borrows the pooled model from all selected chromosomes,
not an autosomes-only pool. When all cores have local support, no pooled fit is attempted.
Automatic alpha depends on eligible bases in the selected chromosomes. Detect once over
the intended set and apply the weights to subsets. A common manual alpha removes that
selection dependence only for alpha, not for any pooled fallback parameters.

Interpreting the preservation boundary
--------------------------------------
An isolated contiguous call longer than 5 kb in a fully eligible default 100 kb window
exceeds 5% and is preserved. Exactly 5 kb is equality and does not trigger preservation.
Masks or clipped windows with smaller denominators can lower this length toward 2.5 kb
at the 50 kb denominator floor. Nearby calls or joined flanks can preserve shorter calls.
A biological event's full width is not its called-base count. Scattered threshold crossings
inside a genuine gain can still remain below the cutoff and receive correction.

Failures and troubleshooting
-----------------------------
Analysis or report-writing failures publish no new normal tables, BEDs, PNGs, weights, or
this README. Completed normal files are staged until analysis and reporting finish.
An error while moving completed files to their final names can leave a partial set of new results.

Use `outliers.error.txt` after a failure. It records the failed stage, settings, original
error and companion paths. Fit failures additionally save plain `outliers.error_histogram.tsv`
with each failing context's coverage counts. Global rows pool all selected chromosomes
and have `NA` coordinates. Second-fit errors include T1, alpha and initial parameters.
After the second pass, `--ambiguous-region-action error` reports every failed merged region
across all selected chromosomes in `outliers.error_audit.tsv`. It includes member calls,
selected windows, eligible/called counts, scores, cutoff, and contributing exact calls.
Actions are `withheld_error`. The table covers failed regions only, not successful decisions.
No additional BAM pass or plots are needed to collect these failures.

Error filenames have the same optional prefix and never reuse success filenames. Retries
replace matching error files. Their attempt identifiers distinguish stale companions.
A successful retry may leave old error files and optional PNGs. Earlier successful normal
files also remain intact after an analysis failure, so check the exit status and current
attempt, rather than treating any pre-existing file as evidence of this attempt's success.
If saving diagnostics also fails, the terminal reports both errors.

An initial-fit failure cannot be repaired just by lowering alpha: alpha is used afterward.
Inspect filters, zero coverage and context counts. Increasing local contexts may help local
support but cannot change the pooled histogram. A smaller alpha can retain more positive
coverage levels for an insufficient or unidentified second fit. T1=2 cannot identify either
free-zero mixture. At T1=3, ZINB rejects the interior dispersion/zero-inflation ridge while
retaining valid constrained boundary fits. The ridge means different model parameters fit
the retained counts equally well but predict different full tails. No threshold is raised
or model substituted silently. ZIP is a different assumption, not a universal recovery mode.

Metadata names and counters
---------------------------
`model_core_size` maps to `--stride`, `minimum_context_span` to `--bin-size`,
`minimum_context_fragment_support` to `--min-context-fragments`, `minimum_mapping_quality`
to `--min-mapq`, and `fragment_length_range` to the accepted fragment length limits.
`detector` records the two-stage algorithm as well as family. Automatic-only
`automatic_probability_denominator` and `tail_probability_multiplier` are `NA` with manual alpha. `eligible_positions` remains
available. Names will be reviewed together with the CLI and consumers.
Zero-weight exclusion counters count tile encounters. Padding-only fragments can make these
counts differ across GC modes, without demonstrating a difference in their counted coverage.
The counter discrepancy is documented for later work.

Genomic extent and coverage reduction
------------------------------------
Called-base percentages describe how much genomic sequence was called. They do not
measure the fraction of coverage removed. Whole-fragment weighting also changes
coverage outside exact calls, which is not captured by the audit's mass columns.
Compare raw and weighted `fcoverage` outputs using identical filters and blacklists.
Start without GC correction, genomic scaling, or length normalization to compare the
coverage directly with the model's mean and cutoff. Then assess the correction with
the profiles or features that matter for your analysis.
"#;
