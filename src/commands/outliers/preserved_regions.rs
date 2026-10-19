//! Group preserved calls for inspection without changing preservation or correction decisions.

use super::{config::OutliersConfig, outliers::CalledRegion, preservation::ChromosomeDecisions};
use crate::shared::interval::Interval;
use anyhow::{Context, Result, ensure};

/// A reporting group of preserved calls, with exact-call bounds and a wider plotting context.
///
/// Existing merged preserved regions remain intact. Their selected window spans connect nearby
/// regions into inspection groups. These bounds locate the evidence, not biological breakpoints.
#[derive(Debug)]
pub(super) struct PreservedRegion {
    pub(super) region_id: String,
    pub(super) chromosome: String,
    /// Smallest span containing all exact member calls, including gaps between them.
    pub(super) interval: Interval<u32>,
    /// Bounds containing the selected preservation windows of all member calls.
    pub(super) window_span: Interval<u32>,
    pub(super) plot_interval: Interval<u32>,
    /// Sorted indices into this chromosome's exact candidates, including calls protected by flanks.
    pub(super) call_indices: Vec<usize>,
    /// Sum of exact called lengths, excluding gaps and flanks.
    pub(super) called_bases: u64,
    pub(super) maximum_window_percentage: f64,
}

/// Collect one chromosome's preserved reporting groups in genomic order.
///
/// Overlapping selected window spans join groups transitively. Merely touching spans remain
/// separate. Plot contexts do not participate in grouping, even when those wider contexts overlap.
pub(super) fn collect_preserved_regions(
    chromosome: &str,
    regions: &[CalledRegion],
    decisions: &ChromosomeDecisions,
    config: &OutliersConfig,
    chromosome_length: u32,
    first_region_number: usize,
) -> Result<Vec<PreservedRegion>> {
    ensure!(
        regions.len() == decisions.calls.len(),
        "preserved-region call decisions do not match candidates"
    );
    let mut seeds = Vec::new();
    for merged in decisions
        .merged_regions
        .iter()
        .filter(|merged| merged.preserved)
    {
        // Keep every member of an already preserved merged region in the same reporting group
        let window_span = decisions.calls[merged.members.clone()]
            .iter()
            .map(|decision| decision.window.interval)
            .reduce(|span, window| span.expand_to_include(window))
            .context("preserved merged region has no member calls")?;
        seeds.push((window_span, merged.members.clone().collect::<Vec<_>>()));
    }
    // Backward fallback windows can start before earlier calls' selected windows, so sort by
    // window position before joining rather than assuming that call order is also window order
    seeds.sort_unstable_by_key(|(span, _)| (span.start(), span.end()));
    let mut groups: Vec<(Interval<u32>, Vec<usize>)> = Vec::new();
    for (window_span, members) in seeds {
        if let Some((previous_span, previous_members)) = groups.last_mut()
            && previous_span.intersects(window_span)
        {
            *previous_span = previous_span.expand_to_include(window_span);
            previous_members.extend(members);
        } else {
            groups.push((window_span, members));
        }
    }
    for (_, members) in &mut groups {
        members.sort_unstable();
    }
    groups.sort_unstable_by_key(|(_, members)| regions[members[0]].interval.start());
    groups
        .into_iter()
        .enumerate()
        .map(|(group_index, (window_span, members))| {
            let interval = regions[members[0]].interval.expand_to_include(
                regions[*members
                    .last()
                    .context("preserved reporting group has no calls")?]
                .interval,
            );
            let plot_interval = plot_context(interval, window_span, config, chromosome_length)?;
            Ok(PreservedRegion {
                region_id: format!("region_{:06}", first_region_number + group_index),
                chromosome: chromosome.to_string(),
                interval,
                window_span,
                plot_interval,
                called_bases: members
                    .iter()
                    .map(|&index| u64::from(regions[index].interval.len()))
                    .sum(),
                maximum_window_percentage: members
                    .iter()
                    .map(|&index| decisions.calls[index].window.score)
                    .fold(0.0, f64::max),
                call_indices: members,
            })
        })
        .collect()
}

/// Center the requested context on the called span, include its evidence windows, and clip bounds.
fn plot_context(
    interval: Interval<u32>,
    window_span: Interval<u32>,
    config: &OutliersConfig,
    chromosome_length: u32,
) -> Result<Interval<u32>> {
    // Bound padding, rather than clipping the requested width before centering. A requested
    // view wider than the chromosome should still reach both ends when the group is off-center
    let maximum_useful_span = u64::from(interval.len()) + 2 * u64::from(chromosome_length);
    let requested_span = (f64::from(interval.len()) * config.preserved_plot_context_multiplier)
        .max(f64::from(config.preserved_plot_min_span))
        .min(maximum_useful_span as f64)
        .ceil() as u64;
    let extra = requested_span.saturating_sub(u64::from(interval.len()));
    let left_padding = (extra / 2).min(u64::from(chromosome_length)) as u32;
    let right_padding = (extra - extra / 2).min(u64::from(chromosome_length)) as u32;
    let context = Interval::new(
        interval.start().saturating_sub(left_padding),
        interval
            .end()
            .saturating_add(right_padding)
            .min(chromosome_length),
    )?;
    Ok(context.expand_to_include(window_span))
}

#[cfg(test)]
mod tests {
    include!("preserved_regions_tests.rs");
}
