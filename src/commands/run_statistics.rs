use crate::{
    commands::counters::BaseCounters,
    shared::gc_tag::{
        MAX_REASONABLE_GC_WEIGHT, MIN_REASONABLE_GC_WEIGHT, gc_failure_action_description,
    },
};
use std::time::Duration;

/// Shared note for tile-based commands whose overlap halos can inflate statistics only.
pub(crate) const TILE_DOUBLE_COUNT_NOTE: &str = "Note: A few reads/fragments may be counted twice in the statistics (only) around the parallelization tiles.";

/// Labels for the common fragment-statistics block.
#[derive(Clone, Copy)]
pub(crate) struct FragmentStatisticsLabels<'a> {
    pub total_reads: &'a str,
    pub accepted_reads: &'a str,
    pub counted_fragments: &'a str,
}

/// Default labels for commands that report ordinary fragment counting.
pub(crate) const DEFAULT_FRAGMENT_STATISTICS_LABELS: FragmentStatisticsLabels<'static> =
    FragmentStatisticsLabels {
        total_reads: "Total observed reads",
        accepted_reads: "Initially accepted reads",
        counted_fragments: "Fragments counted one or more times",
    };

/// Optional GC-related statistics shared by commands that use GC correction or GC tags.
#[derive(Clone, Copy)]
pub(crate) struct GCStatisticsSummary {
    pub neutralize_invalid_gc: bool,
    pub failed_fragments: u64,
    pub missing_tags: Option<u64>,
    pub out_of_range_tags: Option<u64>,
}

/// Common options for fragment-processing command statistics.
#[derive(Clone, Copy)]
pub(crate) struct FragmentRunStatisticsOptions<'a> {
    pub include_section_header: bool,
    pub notes: &'a [&'a str],
    pub labels: FragmentStatisticsLabels<'a>,
    pub blacklist_excluded_fragments: Option<u64>,
    pub gc: Option<GCStatisticsSummary>,
}

/// Print the shared statistics block used by fragment-oriented commands.
///
/// This centralizes the repeated read acceptance, optional blacklist, optional GC,
/// and counted-fragment reporting while leaving command-specific extras to the caller.
pub(crate) fn print_fragment_run_statistics<I, S>(
    base: &BaseCounters,
    elapsed: Duration,
    options: FragmentRunStatisticsOptions<'_>,
    extra_lines: I,
) where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    if options.include_section_header {
        println!();
        println!("Statistics");
        println!("----------");
    }

    for note in options.notes {
        println!("  {}", note);
    }

    let accepted_reads = base.accepted_forward + base.accepted_reverse;
    let accepted_pct = if base.total_reads == 0 {
        0.0
    } else {
        accepted_reads as f64 / base.total_reads as f64 * 100.0
    };

    println!("  {}: {}", options.labels.total_reads, base.total_reads);
    println!(
        "  {}: {} ({:.2}%, forward: {}, reverse: {})",
        options.labels.accepted_reads,
        accepted_reads,
        accepted_pct,
        base.accepted_forward,
        base.accepted_reverse
    );

    if let Some(blacklisted_fragments) = options.blacklist_excluded_fragments {
        println!("  Blacklist-excluded fragments: {}", blacklisted_fragments);
    }

    if let Some(gc) = options.gc {
        let gc_fail_action = gc_failure_action_description(gc.neutralize_invalid_gc);
        println!(
            "  GC correction failures ({}): {}",
            gc_fail_action, gc.failed_fragments
        );

        if let Some(missing_tags) = gc.missing_tags
            && missing_tags > 0
        {
            let missing_action = if gc.neutralize_invalid_gc {
                "counted with weight 1.0 via --neutralize-invalid-gc"
            } else {
                "skipped by default"
            };
            println!(
                "  Warning: fragments missing GC tags: {} ({})",
                missing_tags, missing_action
            );
        }

        if let Some(out_of_range_tags) = gc.out_of_range_tags
            && out_of_range_tags > 0
        {
            println!(
                "  Non-zero GC tag values outside the supported positive range [{:.0e}, {:.0e}] treated as invalid: {}",
                MIN_REASONABLE_GC_WEIGHT, MAX_REASONABLE_GC_WEIGHT, out_of_range_tags
            );
        }
    }

    println!(
        "  {}: {}",
        options.labels.counted_fragments, base.counted_fragments
    );

    for line in extra_lines {
        println!("  {}", line.as_ref());
    }

    if options.include_section_header {
        println!("----------");
    }
    println!("Elapsed time: {:.2?}", elapsed);
}
