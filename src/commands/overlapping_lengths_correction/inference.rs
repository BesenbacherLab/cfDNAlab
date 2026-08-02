//! Bounded rolling construction and application of positional overlap-length correction bins.

use std::{collections::BTreeMap, sync::Arc};

use anyhow::{Context, Result, ensure};

use crate::shared::{fragment::segment_fragment::FragmentWithSegments, interval::Interval};

use super::package::OverlappingLengthsCorrectionPackage;

/// Bin marker used where raw overlap coverage is zero and no average fragment length exists.
const NO_OVERLAP_LENGTH_BIN: u32 = u32::MAX;

/// Difference-array event at a covered segment boundary.
#[derive(Clone, Copy, Debug, Default)]
struct OverlapEvent {
    /// Change in the number of fragments covering the next reference base.
    depth_delta: i64,
    /// Change in the sum of full directional fragment lengths covering the next reference base.
    fragment_length_sum_delta: i64,
}

/// Tile-core lookup bins for positional average overlapping fragment length correction.
///
/// The array stores bin indices instead of `f64` multipliers to keep per-worker memory bounded.
/// Applying the correction reads one compact index per position. Consecutive positions in the same
/// bin reuse the preceding package multiplier lookup.
#[derive(Debug)]
pub(crate) struct PositionalOverlapLengthBins {
    bin_indices: Vec<u32>,
    package: Arc<OverlappingLengthsCorrectionPackage>,
}

impl PositionalOverlapLengthBins {
    /// Multiply finalized tile-core coverage by its positional overlap-length correction.
    ///
    /// Positions without raw overlap coverage retain their existing value. In a consistent
    /// fcoverage run those positions already have zero coverage, but treating the marker as neutral
    /// keeps this helper explicit and prevents an invalid array access.
    pub(crate) fn apply_to_coverage(&self, coverage: &mut [f32]) -> Result<()> {
        ensure!(
            coverage.len() == self.bin_indices.len(),
            "positional overlap-length bins contain {} bases but coverage contains {}",
            self.bin_indices.len(),
            coverage.len()
        );

        let mut previous_bin_index = NO_OVERLAP_LENGTH_BIN;
        let mut previous_weight = 1.0_f64;
        for (position, (&bin_index, coverage_value)) in
            self.bin_indices.iter().zip(coverage.iter_mut()).enumerate()
        {
            if bin_index == NO_OVERLAP_LENGTH_BIN {
                continue;
            }
            if bin_index != previous_bin_index {
                previous_weight = self.package.weight_for_bin_index(bin_index)?;
                previous_bin_index = bin_index;
            }
            let corrected_value = f64::from(*coverage_value) * previous_weight;
            ensure!(
                corrected_value.is_finite() && corrected_value <= f64::from(f32::MAX),
                "overlap-length correction produced invalid coverage at tile position {}: {}",
                position,
                corrected_value
            );
            *coverage_value = corrected_value as f32;
        }
        Ok(())
    }

    #[cfg(test)]
    fn bin_indices(&self) -> &[u32] {
        &self.bin_indices
    }
}

/// Collect positional overlap context while borrowing fragments from the normal fcoverage stream.
///
/// The collector tracks the largest observed `fragment.start()`. A later fragment may begin up to
/// `maximum_fragment_length` bases before that value, so positions earlier than
/// `largest_seen_start - maximum_fragment_length` are complete. Each complete constant-context run
/// fills the corresponding part of the tile-core bin array. Fragments remain owned by the caller
/// and can therefore proceed directly through GC correction and ordinary coverage accumulation.
pub(crate) struct PositionalOverlapLengthCollector {
    /// Loaded correction lookup used to assign average fragment lengths to bins.
    package: Arc<OverlappingLengthsCorrectionPackage>,
    /// Inclusive fragment length bound used to prove that earlier positions are complete.
    maximum_fragment_length: u32,
    /// Reference interval in which fragment overlap events may be collected.
    context: Interval<u32>,
    /// Non-overlapping tile core receiving the final bin indices.
    core: Interval<u32>,
    /// Sparse start and end changes for borrowed fragment segments.
    events: BTreeMap<u32, OverlapEvent>,
    /// Bin index for every tile-core base.
    bin_indices: Vec<u32>,
    /// First genomic position whose overlap context is not yet complete.
    finalized_until: u32,
    /// Raw fragment count at `finalized_until`.
    active_depth: i64,
    /// Sum of full directional fragment lengths at `finalized_until`.
    active_fragment_length_sum: i64,
    /// Largest directional fragment start observed so far.
    largest_seen_start: Option<u32>,
}

impl PositionalOverlapLengthCollector {
    /// Create a rolling collector for one tile core and its fetched context.
    ///
    /// The incoming fragment stream must already enforce `maximum_fragment_length` and originate
    /// from cfDNAlab's coordinate-sorted pairing path. Windowed fcoverage may narrow `context` to
    /// only the requested part of the core. Core positions outside that context retain the neutral
    /// marker because they are not part of the requested output.
    pub(crate) fn new(
        package: Arc<OverlappingLengthsCorrectionPackage>,
        context: Interval<u32>,
        core: Interval<u32>,
    ) -> Result<Self> {
        let core_length = usize::try_from(core.len()).context("tile core length exceeds usize")?;
        let maximum_fragment_length = package.maximum_fragment_length;
        Ok(Self {
            package,
            maximum_fragment_length,
            context,
            core,
            events: BTreeMap::new(),
            bin_indices: vec![NO_OVERLAP_LENGTH_BIN; core_length],
            finalized_until: context.start(),
            active_depth: 0,
            active_fragment_length_sum: 0,
            largest_seen_start: None,
        })
    }

    /// Borrow one raw accepted fragment and update positional overlap context.
    ///
    /// This must be called before any later GC validation can reject the fragment. Segment-aware
    /// fragments update only their counted reference segments, while each segment contributes the
    /// fragment's full directional `forward.pos` to `reverse.reference_end` length.
    pub(crate) fn observe(&mut self, fragment: &FragmentWithSegments) -> Result<()> {
        self.largest_seen_start = Some(
            self.largest_seen_start
                .map_or(fragment.start(), |largest| largest.max(fragment.start())),
        );
        if let Some(segments) = &fragment.segments {
            for segment in segments {
                self.add_segment_events(*segment, fragment.len())?;
            }
        } else {
            self.add_segment_events(fragment.interval, fragment.len())?;
        }
        self.finalize_safe_positions()
    }

    /// Finish the remaining reference context and return the tile-core bin array.
    pub(crate) fn finish(mut self) -> Result<PositionalOverlapLengthBins> {
        self.finalize_until(self.context.end())?;
        Ok(PositionalOverlapLengthBins {
            bin_indices: self.bin_indices,
            package: self.package,
        })
    }

    /// Add clipped start and end events for one counted segment.
    fn add_segment_events(&mut self, segment: Interval<u32>, fragment_length: u32) -> Result<()> {
        let Some(segment) = segment.clip_to(self.context) else {
            return Ok(());
        };
        ensure!(
            segment.start() >= self.finalized_until,
            "fragment beginning at {} arrived after overlap-length position {} was finalized",
            segment.start(),
            self.finalized_until
        );
        let start_event = self.events.entry(segment.start()).or_default();
        start_event.depth_delta += 1;
        start_event.fragment_length_sum_delta += i64::from(fragment_length);
        let end_event = self.events.entry(segment.end()).or_default();
        end_event.depth_delta -= 1;
        end_event.fragment_length_sum_delta -= i64::from(fragment_length);
        Ok(())
    }

    /// Finalize every position that no future accepted fragment can overlap.
    fn finalize_safe_positions(&mut self) -> Result<()> {
        let Some(largest_seen_start) = self.largest_seen_start else {
            return Ok(());
        };
        let safe_until = largest_seen_start
            .saturating_sub(self.maximum_fragment_length)
            .min(self.context.end());
        self.finalize_until(safe_until)
    }

    /// Resolve sparse events and write one bin value for each complete constant-context run.
    fn finalize_until(&mut self, safe_until: u32) -> Result<()> {
        if safe_until <= self.finalized_until {
            return Ok(());
        }
        while self.finalized_until < safe_until {
            if let Some(event) = self.events.remove(&self.finalized_until) {
                self.active_depth += event.depth_delta;
                self.active_fragment_length_sum += event.fragment_length_sum_delta;
            }
            ensure!(
                self.active_depth >= 0 && self.active_fragment_length_sum >= 0,
                "rolling overlap-length context became negative at position {}",
                self.finalized_until
            );

            let run_end = self
                .events
                .first_key_value()
                .map(|(&position, _)| position)
                .unwrap_or(safe_until)
                .min(safe_until);
            ensure!(
                run_end > self.finalized_until,
                "overlap-length event ordering did not advance beyond position {}",
                self.finalized_until
            );
            self.fill_core_bin_indices(self.finalized_until, run_end)?;
            self.finalized_until = run_end;
        }
        Ok(())
    }

    /// Fill the tile-core part of a constant raw depth and fragment length sum run.
    fn fill_core_bin_indices(&mut self, run_start: u32, run_end: u32) -> Result<()> {
        if self.active_depth == 0 || run_end <= self.core.start() || run_start >= self.core.end() {
            return Ok(());
        }
        ensure!(
            self.active_fragment_length_sum > 0,
            "covered overlap-length run {}-{} has a zero fragment length sum",
            run_start,
            run_end
        );
        let average_length = self.active_fragment_length_sum as f64 / self.active_depth as f64;
        let bin_index = self.package.bin_index_for_average_length(average_length)?;
        let bin_index = u32::try_from(bin_index)
            .context("overlap-length model contains more bins than u32 can index")?;
        ensure!(
            bin_index != NO_OVERLAP_LENGTH_BIN,
            "overlap-length model bin index conflicts with the uncovered-position marker"
        );

        let clipped_start = run_start.max(self.core.start());
        let clipped_end = run_end.min(self.core.end());
        let local_start = (clipped_start - self.core.start()) as usize;
        let local_end = (clipped_end - self.core.start()) as usize;
        self.bin_indices[local_start..local_end].fill(bin_index);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    include!("inference_tests.rs");
}
