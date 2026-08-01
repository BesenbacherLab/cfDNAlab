//! Bounded rolling inference of scalar overlap-length weights for the fcoverage fragment stream.

use std::{
    collections::{BTreeMap, VecDeque},
    sync::Arc,
};

use anyhow::{Context, Result, ensure};

use crate::shared::{fragment::segment_fragment::FragmentWithSegments, interval::Interval};

use super::package::OverlappingLengthsCorrectionPackage;

/// Number of finalized reference bases stored in each rolling prefix chunk.
///
/// A fixed chunk size keeps allocation and lookup straightforward while allowing old coordinate
/// ranges to be discarded without moving the retained prefix values.
const PREFIX_CHUNK_BASES: usize = 65_536;

/// Difference-array event at a covered segment boundary.
#[derive(Clone, Copy, Debug, Default)]
struct OverlapEvent {
    /// Change in the number of fragments covering the next reference base.
    depth_delta: i64,
    /// Change in the sum of full directional fragment lengths covering the next reference base.
    fragment_length_sum_delta: i64,
}

/// Prefix sums for a contiguous chunk of finalized genomic positions.
///
/// Both vectors contain an initial value at `start`, followed by one value for every finalized base
/// in the chunk. Seeding a new chunk with the preceding cumulative values makes either endpoint of
/// a segment queryable without retaining earlier chunks.
#[derive(Debug)]
struct PrefixChunk {
    /// Genomic coordinate represented by prefix offset zero.
    start: u32,
    /// Cumulative sum of positional correction weights.
    weight_prefix: Vec<f64>,
    /// Cumulative count of non-blacklisted covered positions.
    eligible_prefix: Vec<u64>,
}

impl PrefixChunk {
    /// Start a chunk with cumulative values inherited from the preceding chunk.
    fn new(start: u32, preceding_weight: f64, preceding_eligible: u64) -> Self {
        Self {
            start,
            weight_prefix: vec![preceding_weight],
            eligible_prefix: vec![preceding_eligible],
        }
    }

    /// Number of genomic positions finalized into this chunk.
    fn position_count(&self) -> usize {
        self.weight_prefix.len() - 1
    }

    /// Exclusive genomic end of the positions stored in this chunk.
    fn end(&self) -> u32 {
        self.start + self.position_count() as u32
    }

    /// Return whether this chunk can answer a prefix query at `coordinate`.
    ///
    /// The exclusive position end is included because prefix arrays have one more value than the
    /// number of represented bases.
    fn contains_prefix_coordinate(&self, coordinate: u32) -> bool {
        self.start <= coordinate && coordinate <= self.end()
    }
}

/// Add average-overlapping-length weights while retaining only bounded positional state.
///
/// The iterator tracks the largest `fragment.start()` returned by the inner fragment iterator. A
/// later fragment may start up to `maximum_fragment_length` bases before that value, so positions
/// earlier than `largest_seen_start - maximum_fragment_length` are safe to finalize. Finalized
/// correction prefixes are kept in 64 KiB chunks and old chunks are discarded as soon as no queued
/// fragment can refer to them.
pub(crate) struct OverlappingLengthWeightIterator<I> {
    /// Normal fragment iterator whose exact return order must be preserved.
    inner: I,
    /// Loaded correction lookup, or `None` for direct pass-through mode.
    package: Option<Arc<OverlappingLengthsCorrectionPackage>>,
    /// Inclusive fragment length bound used to prove that earlier positions are final.
    maximum_fragment_length: u32,
    /// Reference interval for which overlap context may be accumulated.
    context: Interval<u32>,
    /// Sorted chromosome blacklist intervals used to exclude positional weights.
    blacklist: Vec<Interval<u64>>,
    /// First blacklist interval that may overlap `finalized_until`.
    blacklist_index: usize,
    /// Sparse start and end changes for ingested fragment segments.
    events: BTreeMap<u32, OverlapEvent>,
    /// Ingested fragments in the exact order returned by `inner`.
    pending: VecDeque<FragmentWithSegments>,
    /// Weighted fragments that are safe for the caller to consume.
    ready: VecDeque<FragmentWithSegments>,
    /// Rolling chunks containing finalized positional prefix sums.
    prefix_chunks: VecDeque<PrefixChunk>,
    /// Number of finalized genomic positions stored in each prefix chunk.
    ///
    /// Production construction always uses `PREFIX_CHUNK_BASES`. Keeping the value on the iterator
    /// allows unit tests to prove that changing only the storage partition leaves results
    /// unchanged.
    prefix_chunk_bases: usize,
    /// First genomic position not yet finalized.
    finalized_until: u32,
    /// Raw fragment count covering `finalized_until` during the left-to-right scan.
    active_depth: i64,
    /// Sum of full directional fragment lengths covering `finalized_until`.
    active_fragment_length_sum: i64,
    /// Cumulative correction-weight sum at `finalized_until`.
    cumulative_weight: f64,
    /// Cumulative eligible-base count at `finalized_until`.
    cumulative_eligible: u64,
    /// Largest directional fragment start returned by `inner` so far.
    largest_seen_start: Option<u32>,
    /// Whether the inner iterator has ended or returned an error.
    reached_end: bool,
}

impl<I> OverlappingLengthWeightIterator<I>
where
    I: Iterator<Item = Result<FragmentWithSegments>>,
{
    /// Wrap a fragment iterator with bounded overlap-context calculation.
    ///
    /// `context` must include every base needed to weight fragments that can contribute to the tile
    /// core. fcoverage supplies a two-maximum-fragment-length tile halo so the inner span contains
    /// both contributing fragments and the neighboring fragments defining their overlap context.
    /// The inner iterator must yield only fragments no longer than `maximum_fragment_length` and
    /// must originate from cfDNAlab's coordinate-sorted pairing adaptor. Those conditions make the
    /// largest-seen-start finalization boundary valid.
    pub(crate) fn new(
        inner: I,
        package: Option<Arc<OverlappingLengthsCorrectionPackage>>,
        maximum_fragment_length: u32,
        context: Interval<u32>,
        blacklist: &[Interval<u64>],
    ) -> Self {
        // A pass-through iterator never reads or retains positional masking state
        let (blacklist, blacklist_index) = if package.is_some() {
            (
                blacklist.to_vec(),
                blacklist.partition_point(|interval| interval.end() <= u64::from(context.start())),
            )
        } else {
            (Vec::new(), 0)
        };
        Self {
            inner,
            package,
            maximum_fragment_length,
            context,
            blacklist,
            blacklist_index,
            events: BTreeMap::new(),
            pending: VecDeque::new(),
            ready: VecDeque::new(),
            prefix_chunks: VecDeque::new(),
            prefix_chunk_bases: PREFIX_CHUNK_BASES,
            finalized_until: context.start(),
            active_depth: 0,
            active_fragment_length_sum: 0,
            cumulative_weight: 0.0,
            cumulative_eligible: 0,
            largest_seen_start: None,
            reached_end: false,
        }
    }

    /// Replace the production chunk size for storage-partition equivalence tests.
    ///
    /// This method is unavailable in production builds. A positive chunk size is required because
    /// a zero-sized chunk could never accept a finalized genomic position.
    #[cfg(test)]
    fn with_prefix_chunk_bases(mut self, prefix_chunk_bases: usize) -> Self {
        assert!(prefix_chunk_bases > 0, "prefix chunk size must be positive");
        self.prefix_chunk_bases = prefix_chunk_bases;
        self
    }

    /// Recover the normal fragment iterator after this adaptor has been drained.
    ///
    /// fcoverage uses this to read the pairing iterator's local counters without requiring the
    /// overlap adaptor to know anything about counter implementations.
    pub(crate) fn into_inner(self) -> I {
        self.inner
    }

    /// Add a returned fragment to overlap events and the pending-fragment queue.
    ///
    /// Segment-aware fragments update only their counted reference segments. Every updated base
    /// nevertheless receives the fragment's full directional `forward.pos` to
    /// `reverse.reference_end` length.
    fn ingest(&mut self, fragment: FragmentWithSegments) -> Result<()> {
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
        self.pending.push_back(fragment);
        Ok(())
    }

    /// Add clipped start and end events for a single counted segment.
    ///
    /// A segment starting before `finalized_until` would invalidate the bounded-stream proof and is
    /// treated as an ordering error rather than silently changing an already assigned weight.
    fn add_segment_events(&mut self, segment: Interval<u32>, fragment_length: u32) -> Result<()> {
        let Some(segment) = segment.clip_to(self.context) else {
            return Ok(());
        };
        ensure!(
            segment.start() >= self.finalized_until,
            "fragment beginning at {} arrived after overlapping-length position {} was finalized",
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
    ///
    /// Returned fragment starts need not be sorted because pairs are returned when their second
    /// mate is consumed. Nevertheless, after observing a largest start `S`, a future accepted
    /// fragment of maximum length `M` cannot begin before `S - M`. Positions before that coordinate
    /// are therefore complete.
    fn finalize_safe_positions(&mut self) -> Result<()> {
        let Some(largest_seen_start) = self.largest_seen_start else {
            return Ok(());
        };
        let safe_until = largest_seen_start
            .saturating_sub(self.maximum_fragment_length)
            .min(self.context.end());
        self.finalize_until(safe_until)
    }

    /// Resolve overlap events and append positional correction prefixes up to `safe_until`.
    ///
    /// `safe_until` is exclusive. At each covered, non-blacklisted base, raw depth and the full
    /// fragment length sum define an average length used for package lookup. Uncovered or
    /// blacklisted bases contribute neither weight nor eligible-base count.
    fn finalize_until(&mut self, safe_until: u32) -> Result<()> {
        if safe_until <= self.finalized_until {
            return Ok(());
        }
        let package = self
            .package
            .as_ref()
            .context("overlapping-length package missing while finalizing correction positions")?
            .clone();
        // Resolve sparse difference events in genomic order without a chromosome-sized array
        while self.finalized_until < safe_until {
            let position = self.finalized_until;
            if let Some(event) = self.events.remove(&position) {
                self.active_depth += event.depth_delta;
                self.active_fragment_length_sum += event.fragment_length_sum_delta;
            }
            ensure!(
                self.active_depth >= 0 && self.active_fragment_length_sum >= 0,
                "rolling overlapping fragment length prefix became negative at position {}",
                position
            );
            // Advance the sorted blacklist cursor once intervals are entirely behind the scan
            while self.blacklist_index < self.blacklist.len()
                && self.blacklist[self.blacklist_index].end() <= u64::from(position)
            {
                self.blacklist_index += 1;
            }
            let blacklisted = self.blacklist_index < self.blacklist.len()
                && self.blacklist[self.blacklist_index].start() <= u64::from(position)
                && u64::from(position) < self.blacklist[self.blacklist_index].end();
            // Lookup uses raw overlap context even when fcoverage later applies GC weights
            let (weight, eligible) = if self.active_depth > 0 && !blacklisted {
                let average_length =
                    self.active_fragment_length_sum as f64 / self.active_depth as f64;
                (package.weight_for_average_length(average_length)?, 1_u64)
            } else {
                (0.0, 0_u64)
            };
            self.append_prefix_position(position, weight, eligible);
            self.finalized_until += 1;
            // Chunk boundaries are natural opportunities to release fragments and reclaim memory
            if self
                .prefix_chunks
                .back()
                .is_some_and(|chunk| chunk.position_count() == self.prefix_chunk_bases)
            {
                self.release_finalized_fragments()?;
                self.drop_expired_chunks();
            }
        }
        self.release_finalized_fragments()?;
        Ok(())
    }

    /// Append one finalized position to the rolling cumulative prefixes.
    ///
    /// A new chunk inherits the preceding cumulative totals. Segment sums can therefore subtract
    /// endpoints from different retained chunks without any adjustment.
    fn append_prefix_position(&mut self, position: u32, weight: f64, eligible: u64) {
        let needs_chunk = self
            .prefix_chunks
            .back()
            .is_none_or(|chunk| chunk.position_count() == self.prefix_chunk_bases);
        if needs_chunk {
            self.prefix_chunks.push_back(PrefixChunk::new(
                position,
                self.cumulative_weight,
                self.cumulative_eligible,
            ));
        }
        self.cumulative_weight += weight;
        self.cumulative_eligible += eligible;
        let chunk = self
            .prefix_chunks
            .back_mut()
            .expect("prefix chunk was inserted above");
        chunk.weight_prefix.push(self.cumulative_weight);
        chunk.eligible_prefix.push(self.cumulative_eligible);
    }

    /// Assign weights to the consecutive ready fragments at the front of the FIFO queue.
    ///
    /// Only the oldest queued fragment is inspected. If it is not ready, later fragments remain
    /// queued even when their shorter spans are already finalized. This makes the adaptor's output
    /// order exactly equal to the inner fragment iterator's output order.
    fn release_finalized_fragments(&mut self) -> Result<()> {
        while let Some(fragment) = self.pending.front() {
            if fragment.end() > self.finalized_until {
                break;
            }
            let mut fragment = self
                .pending
                .pop_front()
                .expect("pending front existed above");
            fragment.overlap_length_weight = self.weight_for_fragment(&fragment)?;
            self.ready.push_back(fragment);
        }
        Ok(())
    }

    /// Calculate a fragment's scalar correction over its original counted span.
    ///
    /// The scalar is the base-pair-weighted mean of positional lookup weights. Blacklisted and
    /// uncovered positions are absent from both numerator and denominator. A fragment with no
    /// eligible bases receives neutral weight one.
    fn weight_for_fragment(&self, fragment: &FragmentWithSegments) -> Result<f64> {
        let mut weight_sum = 0.0;
        let mut eligible_bases = 0_u64;
        if let Some(segments) = &fragment.segments {
            for segment in segments {
                let (segment_weight, segment_eligible) = self.prefix_sum_for_segment(*segment)?;
                weight_sum += segment_weight;
                eligible_bases += segment_eligible;
            }
        } else {
            let (segment_weight, segment_eligible) =
                self.prefix_sum_for_segment(fragment.interval)?;
            weight_sum += segment_weight;
            eligible_bases += segment_eligible;
        }
        if eligible_bases == 0 {
            return Ok(1.0);
        }
        let weight = weight_sum / eligible_bases as f64;
        ensure!(
            weight.is_finite() && weight > 0.0,
            "fragment overlapping-length weight must be finite and positive"
        );
        Ok(weight)
    }

    /// Return correction-weight and eligible-base sums for a clipped counted segment.
    fn prefix_sum_for_segment(&self, segment: Interval<u32>) -> Result<(f64, u64)> {
        let Some(segment) = segment.clip_to(self.context) else {
            return Ok((0.0, 0));
        };
        ensure!(
            segment.end() <= self.finalized_until,
            "fragment overlap weight requested before its complete span was finalized"
        );
        let (end_weight, end_eligible) = self.prefix_at(segment.end())?;
        let (start_weight, start_eligible) = self.prefix_at(segment.start())?;
        Ok((end_weight - start_weight, end_eligible - start_eligible))
    }

    /// Read cumulative correction totals at a finalized genomic coordinate.
    ///
    /// Chunks are searched newest-first because adjacent chunks both contain the shared boundary
    /// prefix coordinate.
    fn prefix_at(&self, coordinate: u32) -> Result<(f64, u64)> {
        for chunk in self.prefix_chunks.iter().rev() {
            if chunk.contains_prefix_coordinate(coordinate) {
                let offset = (coordinate - chunk.start) as usize;
                return Ok((chunk.weight_prefix[offset], chunk.eligible_prefix[offset]));
            }
        }
        anyhow::bail!(
            "rolling overlapping-length prefix for coordinate {} is no longer retained",
            coordinate
        )
    }

    /// Discard prefix chunks that no still-pending fragment can reference.
    ///
    /// Fragment starts are not necessarily ordered inside the FIFO queue, so the retention boundary
    /// is the smallest start across all queued fragments rather than the queue front's start. This
    /// scan runs only at 64 KiB chunk boundaries, not for every fragment. At least one chunk is
    /// retained so a new chunk can inherit the current cumulative totals.
    fn drop_expired_chunks(&mut self) {
        let retain_from = self
            .pending
            .iter()
            .map(FragmentWithSegments::start)
            .min()
            .unwrap_or(self.finalized_until);
        while self.prefix_chunks.len() > 1
            && self
                .prefix_chunks
                .front()
                .is_some_and(|chunk| chunk.end() <= retain_from)
        {
            self.prefix_chunks.pop_front();
        }
    }

    /// Finalize the remaining context and weight every pending fragment at end of input.
    ///
    /// End of the fragment stream proves that no unseen fragment remains, so the complete bounded
    /// context can be resolved immediately. Remaining fragments are weighted and moved to `ready`
    /// strictly from the front of the FIFO queue.
    fn flush_end_of_stream(&mut self) -> Result<()> {
        self.finalize_until(self.context.end())?;
        while let Some(mut fragment) = self.pending.pop_front() {
            fragment.overlap_length_weight = self.weight_for_fragment(&fragment)?;
            self.ready.push_back(fragment);
        }
        Ok(())
    }
}

impl<I> Iterator for OverlappingLengthWeightIterator<I>
where
    I: Iterator<Item = Result<FragmentWithSegments>>,
{
    type Item = Result<FragmentWithSegments>;

    /// Return the next fragment as soon as its complete correction span is finalized.
    ///
    /// With no package, this is a direct pass-through. With correction enabled, each call first
    /// drains already weighted fragments, then consumes normal fragments until the largest returned
    /// start makes at least one consecutive FIFO-front fragment safe. EOF finalizes and drains all
    /// remaining fragments.
    fn next(&mut self) -> Option<Self::Item> {
        if self.package.is_none() {
            return self.inner.next();
        }
        loop {
            // Preserve streaming behavior by returning ready work before reading more BAM records
            if let Some(fragment) = self.ready.pop_front() {
                return Some(Ok(fragment));
            }
            if self.reached_end {
                return None;
            }
            // The normal iterator may return fragment starts in a locally decreasing order
            match self.inner.next() {
                Some(Ok(fragment)) => {
                    if let Err(error) = self
                        .ingest(fragment)
                        .and_then(|_| self.finalize_safe_positions())
                    {
                        self.reached_end = true;
                        return Some(Err(error));
                    }
                }
                Some(Err(error)) => {
                    self.reached_end = true;
                    return Some(Err(error));
                }
                None => {
                    self.reached_end = true;
                    // EOF proves that no future fragment can change any retained position
                    if let Err(error) = self.flush_end_of_stream() {
                        return Some(Err(error));
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    include!("inference_tests.rs");
}
