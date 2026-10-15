//! Physical call-start windows, the isolated-boundary fallback, and shared merged-region actions.

use super::{
    config::{AmbiguousRegionAction, OutliersConfig},
    outliers::CalledRegion,
};
use crate::shared::interval::{Interval, TouchingMergePolicy, merge_intervals};
use anyhow::{Result, ensure};
use std::{cmp::Reverse, collections::BinaryHeap, ops::Range};

/// Selected evidence for a call. Its numerator always uses known exact called intervals.
#[derive(Debug, Clone, Copy)]
pub(super) struct WindowEvidence {
    /// Physical chromosome span of the selected forward or backward preservation window.
    pub(super) interval: Interval<u32>,
    /// Unblacklisted bases in this window, including positions with zero coverage.
    pub(super) eligible_positions: u64,
    /// Exact T2-called bases inside this window, counted before flanking or preservation.
    pub(super) called_positions: u64,
    /// Percentage of called bases using the configured minimum eligible-base denominator.
    pub(super) score: f64,
    /// Whether this call's decision uses its supported backward fallback window.
    pub(super) backward: bool,
}

/// Per-call evidence and index of the flanked merged region sharing its final action.
#[derive(Debug)]
pub(super) struct CallDecision {
    pub(super) window: WindowEvidence,
    /// Whether this call's selected score strictly exceeds the preservation cutoff.
    pub(super) flagged: bool,
    /// Index into this chromosome's merged flanked regions, which determine the final action.
    pub(super) merged_region_index: usize,
}

/// All members retain individual weights but share preservation or correction.
#[derive(Debug)]
pub(super) struct MergedRegion {
    /// Union span after flanking all member calls and merging overlapping/touching flanks.
    pub(super) interval: Interval<u32>,
    /// Indices of the exact calls in this region, using their original chromosome order.
    pub(super) members: Range<usize>,
    /// Whether any member call requests preservation, regardless of the configured action.
    pub(super) flagged: bool,
    /// Whether the selected action actually withholds correction for all member calls.
    pub(super) preserved: bool,
    /// Largest selected window score among member calls, without recomputing it over flanks.
    pub(super) maximum_score: f64,
}

/// Preservation decisions for one chromosome, with calls in the same order as the input.
#[derive(Debug)]
pub(super) struct ChromosomeDecisions {
    pub(super) calls: Vec<CallDecision>,
    pub(super) merged_regions: Vec<MergedRegion>,
}

/// Prefix lengths over a union of intervals, avoiding chromosome-sized base arrays.
pub(super) struct IntervalLengths {
    intervals: Vec<Interval<u64>>,
    prefix_lengths: Vec<u64>,
}

impl IntervalLengths {
    /// Merge the supplied spans so overlapping bases contribute only once to window counts.
    pub(super) fn new(intervals: Vec<Interval<u64>>) -> Self {
        let intervals = merge_intervals(intervals, TouchingMergePolicy::MergeTouching);
        let mut prefix_lengths = Vec::with_capacity(intervals.len() + 1);
        prefix_lengths.push(0);
        for interval in &intervals {
            prefix_lengths.push(prefix_lengths.last().copied().unwrap_or(0) + interval.len());
        }
        Self {
            intervals,
            prefix_lengths,
        }
    }

    /// Count interval-union bases strictly before an absolute chromosome position.
    fn before(&self, position: u64) -> u64 {
        let complete = self
            .intervals
            .partition_point(|interval| interval.end() <= position);
        let partial = self
            .intervals
            .get(complete)
            .map_or(0, |interval| position.saturating_sub(interval.start()));
        self.prefix_lengths[complete] + partial
    }

    /// Count interval-union bases in a physical `[start, end)` window on this chromosome.
    pub(super) fn within(&self, interval: Interval<u32>) -> u64 {
        self.before(interval.end() as u64) - self.before(interval.start() as u64)
    }
}

/// Calculate window evidence from exact calls and blacklisted spans, without reading coverage.
fn window_evidence(
    interval: Interval<u32>,
    called: &IntervalLengths,
    masked: &IntervalLengths,
    minimum_positions: u32,
    backward: bool,
) -> WindowEvidence {
    // Blacklists reduce support within the physical window, without extending its span
    // The denominator floor limits scores at chromosome ends and beside large masked blocks
    let eligible_positions = interval.len() as u64 - masked.within(interval);
    let called_positions = called.within(interval);
    WindowEvidence {
        interval,
        eligible_positions,
        called_positions,
        score: 100.0 * called_positions as f64
            / eligible_positions.max(minimum_positions as u64) as f64,
        backward,
    }
}

/// Evaluate all forward windows, use a sufficiently supported backward fallback when needed,
/// then flank every candidate before assigning the shared merged-region action.
///
/// Input calls are the chromosome's sorted, joined exact intervals. The default 100 kb windows
/// measure their local genomic extent. They do not replace the broader context used to fit T2.
pub(super) fn evaluate_regions(
    config: &OutliersConfig,
    regions: &[CalledRegion],
    blacklist: &[Interval<u64>],
    chromosome_length: u32,
) -> Result<ChromosomeDecisions> {
    let called = IntervalLengths::new(
        regions
            .iter()
            .map(|region| region.interval.try_to_u64())
            .collect::<crate::Result<Vec<_>>>()?,
    );
    let masked = IntervalLengths::new(blacklist.to_vec());
    // Start one physical window at each exact call, using the complete candidate set for counts
    let windows = regions
        .iter()
        .map(|region| {
            let interval = Interval::new(
                region.interval.start(),
                region
                    .interval
                    .start()
                    .saturating_add(config.protection_window_size)
                    .min(chromosome_length),
            )?;
            Ok(window_evidence(
                interval,
                &called,
                &masked,
                config.min_protection_positions,
                false,
            ))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut calls = Vec::with_capacity(regions.len());
    // Retain the highest-scoring forward window overlapping the current call, including windows
    // started by earlier calls. Later calls cannot start inside an already joined exact call
    let mut active = BinaryHeap::new();
    // Furthest end of an earlier forward window with enough eligible support, independent of score
    let mut supported_window_end = 0_u32;
    for (call_index, region) in regions.iter().enumerate() {
        let own = windows[call_index];
        // Nonnegative finite f64 scores have the same ordering as their bit representations.
        // Reverse index breaks score ties in favor of the earliest evaluated window
        active.push((own.score.to_bits(), Reverse(call_index)));
        while let Some(&(_, Reverse(index))) = active.peek() {
            if windows[index].interval.end() > region.interval.start() {
                break;
            }
            active.pop();
        }
        let &(_, Reverse(best_index)) = active
            .peek()
            .expect("the current call supplies an overlapping window");
        let mut selected = windows[best_index];
        let has_supported_earlier_window = supported_window_end > region.interval.start();
        if own.eligible_positions < config.min_protection_positions as u64
            && !has_supported_earlier_window
        {
            // A short forward window can still receive enough support by looking left from the
            // call's end. Known exact calls supply the numerator without another coverage pass
            let backward_interval = Interval::new(
                region
                    .interval
                    .end()
                    .saturating_sub(config.protection_window_size),
                region.interval.end(),
            )?;
            let backward = window_evidence(
                backward_interval,
                &called,
                &masked,
                config.min_protection_positions,
                true,
            );
            if backward.eligible_positions >= config.min_protection_positions as u64 {
                // A supported fallback replaces undersupported forward evidence for this call.
                // It does not create another ordinary window or recursively change earlier calls
                selected = backward;
            }
        }
        if own.eligible_positions >= config.min_protection_positions as u64 {
            supported_window_end = supported_window_end.max(own.interval.end());
        }
        calls.push(CallDecision {
            window: selected,
            flagged: selected.score > config.max_outlier_percentage,
            merged_region_index: 0,
        });
    }

    let mut merged_regions: Vec<MergedRegion> = Vec::new();
    // Include every candidate when grouping flanks. A flagged member can protect neighbors
    // through a chain of touching flanks, while each exact call keeps its own proposed weight
    for (call_index, region) in regions.iter().enumerate() {
        let flanked = Interval::new(
            region.interval.start().saturating_sub(
                config
                    .blacklist_flank
                    .unwrap_or(config.fragment_lengths.max_fragment_length),
            ),
            region
                .interval
                .end()
                .saturating_add(
                    config
                        .blacklist_flank
                        .unwrap_or(config.fragment_lengths.max_fragment_length),
                )
                .min(chromosome_length),
        )?;
        if let Some(previous) = merged_regions.last_mut()
            && flanked.start() <= previous.interval.end()
        {
            previous.interval = previous.interval.expand_to_include(flanked);
            previous.members.end = call_index + 1;
            previous.flagged |= calls[call_index].flagged;
            previous.maximum_score = previous.maximum_score.max(calls[call_index].window.score);
        } else {
            merged_regions.push(MergedRegion {
                interval: flanked,
                members: call_index..call_index + 1,
                flagged: calls[call_index].flagged,
                preserved: false,
                maximum_score: calls[call_index].window.score,
            });
        }
        calls[call_index].merged_region_index = merged_regions.len() - 1;
    }
    for merged in &mut merged_regions {
        ensure!(
            !merged.flagged || config.ambiguous_region_action != AmbiguousRegionAction::Error,
            "merged outlier region {:?} exceeds max_outlier_percentage={}. ambiguous_region_action=error prevents output publication",
            merged.interval.as_tuple(),
            config.max_outlier_percentage
        );
        merged.preserved =
            merged.flagged && config.ambiguous_region_action == AmbiguousRegionAction::Preserve;
    }
    Ok(ChromosomeDecisions {
        calls,
        merged_regions,
    })
}

#[cfg(test)]
mod tests {
    include!("preservation_tests.rs");
}
