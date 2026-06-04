use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, Record};
use smallvec::SmallVec;

use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::gc_tag::{GCTagValue, combine_gc_tag_values, read_gc_tag_from_record};
use crate::shared::indel_mode::IndelMode;
use crate::shared::interval::{Interval, TouchingMergePolicy, merge_sorted_intervals};

/// Represents a fragment together with the reference segments that are safe for k-mer analysis.
#[derive(Debug, Clone)]
pub(crate) struct FragmentWithKmerSegments {
    pub(crate) interval: Interval<u32>,
    pub(crate) segments: SmallVec<[Interval<u32>; 12]>,
    pub(crate) gc_tag: GCTagValue,
}

impl FragmentWithKmerSegments {
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    #[allow(dead_code)]
    pub(crate) fn end(&self) -> u32 {
        self.interval.end()
    }

    #[inline]
    pub(crate) fn len(&self) -> u32 {
        self.interval.len()
    }
}

/// Captures per-read segment data together with indel information that influences k-mer counting.
#[derive(Debug, Clone)]
pub(crate) struct KmerSegmentedReadInfo {
    pub(crate) tid: i32,
    pub(crate) interval: Interval<u32>,
    pub(crate) is_reverse: bool,
    pub(crate) has_insertion: bool,
    pub(crate) has_deletion: bool,
    pub(crate) leading_insertion: bool,
    pub(crate) trailing_insertion: bool,
    pub(crate) ref_mapped_segments: Vec<(u32, u32)>,
    pub(crate) gc_tag: GCTagValue,
}

impl KmerSegmentedReadInfo {
    /// Return the read's inclusive start on the reference.
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Return the read's exclusive end on the reference.
    #[inline]
    pub(crate) fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Return the read's aligned reference span `[pos, end)`.
    #[inline]
    pub(crate) fn aligned_interval(&self) -> Interval<u32> {
        self.interval
    }

    #[inline]
    pub(crate) fn has_indel(&self) -> bool {
        self.has_insertion || self.has_deletion
    }

    #[inline]
    pub(crate) fn absolute_segments(&self) -> Vec<Interval<u32>> {
        if self.ref_mapped_segments.is_empty() {
            vec![self.aligned_interval()]
        } else {
            self.ref_mapped_segments
                .iter()
                .filter_map(|(off, len)| {
                    let start = self.start().saturating_add(*off);
                    let end = start.saturating_add(*len);
                    Interval::new(start, end).ok()
                })
                .collect()
        }
    }
}

impl KmerSegmentedReadInfo {
    /// Build read metadata, optionally collecting reference segments for indel-aware counting.
    pub(crate) fn from_record(
        r: &Record,
        capture_segments: bool,
        gc_tag: Option<&[u8]>,
    ) -> Result<Self> {
        // First pass: gather flags that drive pairing decisions and mate-gap handling.
        let mut has_insertion = false;
        let mut has_deletion = false;
        let mut leading_insertion = false;
        let mut trailing_insertion = false;
        let mut saw_edge_op = false;

        for c in r.cigar().iter() {
            // Any operation that consumes read or reference sequence tells us something about the
            // outermost aligned bases. `=` and `X` behave like `M` here, and deletions/insertions
            // signal gaps we should respect when merging mate gaps later on.
            let affects_edges = matches!(
                *c,
                Cigar::Match(_)
                    | Cigar::Equal(_)
                    | Cigar::Diff(_)
                    | Cigar::Del(_)
                    | Cigar::RefSkip(_)
                    | Cigar::Ins(_)
            );

            if affects_edges {
                if !saw_edge_op {
                    leading_insertion = matches!(*c, Cigar::Ins(_));
                    saw_edge_op = true;
                }
                trailing_insertion = matches!(*c, Cigar::Ins(_));
            }

            match *c {
                Cigar::Del(_) | Cigar::RefSkip(_) => {
                    has_deletion = true;
                }
                Cigar::Ins(_) => {
                    has_insertion = true;
                }
                _ => {}
            }
        }

        // Second pass: only build explicit segments when the caller asked for indel adjustments.
        let ref_mapped_segments = if capture_segments {
            let mut ref_off: u32 = 0;
            let mut seg_start: Option<u32> = None;
            let mut segments: Vec<(u32, u32)> = Vec::new();

            for c in r.cigar().iter() {
                match *c {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        if seg_start.is_none() {
                            seg_start = Some(ref_off);
                        }
                        ref_off = ref_off.saturating_add(len);
                    }
                    Cigar::Del(len) | Cigar::RefSkip(len) => {
                        if let Some(start) = seg_start.take() {
                            let seg_len = ref_off.saturating_sub(start);
                            if seg_len > 0 {
                                segments.push((start, seg_len));
                            }
                        }
                        ref_off = ref_off.saturating_add(len);
                    }
                    Cigar::Ins(_) => {
                        if let Some(start) = seg_start.take() {
                            let seg_len = ref_off.saturating_sub(start);
                            if seg_len > 0 {
                                segments.push((start, seg_len));
                            }
                        }
                    }
                    Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
            }

            if let Some(start) = seg_start {
                let seg_len = ref_off.saturating_sub(start);
                if seg_len > 0 {
                    segments.push((start, seg_len));
                }
            }

            segments
        } else {
            Vec::new()
        };

        let gc_tag_value = gc_tag
            .map(|tag| read_gc_tag_from_record(r, tag))
            .unwrap_or_default();

        Ok(KmerSegmentedReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            has_insertion,
            has_deletion,
            leading_insertion,
            trailing_insertion,
            ref_mapped_segments,
            gc_tag: gc_tag_value,
        })
    }
}

impl PairOrientable for KmerSegmentedReadInfo {
    #[inline]
    fn tid(&self) -> i32 {
        self.tid
    }

    #[inline]
    fn is_reverse(&self) -> bool {
        self.is_reverse
    }

    #[inline]
    fn pos(&self) -> u32 {
        self.start()
    }
}

/// Quick path for fragments that should be treated as a flat span (no indels or ignoring them):
/// trim offsets, optionally bridge the mate gap, and return the resulting segments without running
/// the full indel bookkeeping.
fn collect_flat_fragment(
    forward: &KmerSegmentedReadInfo,
    reverse: &KmerSegmentedReadInfo,
    include_inter_mate_gap: bool,
    end_offset: u32,
    span_start: u32,
    span_end: u32,
    gc_tag: GCTagValue,
) -> Option<FragmentWithKmerSegments> {
    // Trim fixed offsets from both ends so k-mer contexts avoid edge artifacts.
    // `end_offset` is expected to be small (defaults to 0), so most spans dwarf it.
    // When a fragment is shorter than 2 * end_offset, trimming would invert the span.
    // We thus guard by collapsing `trim_end` to `span_start` in that case so the
    // subsequent check returns None.
    let trim_start = span_start.saturating_add(end_offset);
    let trim_end = if span_end > end_offset {
        span_end - end_offset
    } else {
        span_start
    };
    let trim_window = Interval::new(trim_start, trim_end).ok()?;

    let forward_segment = forward.aligned_interval().clip_to(trim_window);
    let reverse_segment = reverse.aligned_interval().clip_to(trim_window);

    let mates_overlap_or_touch = forward.end() >= reverse.start();
    let mut segments: SmallVec<[Interval<u32>; 12]> = SmallVec::new();

    if include_inter_mate_gap || mates_overlap_or_touch {
        segments.push(trim_window);
    } else {
        if let Some(segment) = forward_segment {
            segments.push(segment);
        }

        if let Some(segment) = reverse_segment {
            if let Some(last) = segments.last_mut() {
                if last.end() >= segment.start() {
                    if segment.end() > last.end() {
                        *last = last.expand_to_include(segment);
                    }
                } else {
                    segments.push(segment);
                }
            } else {
                segments.push(segment);
            }
        }
    }

    if segments.is_empty() {
        return None;
    }

    Some(FragmentWithKmerSegments {
        interval: Interval::new(span_start, span_end).ok()?,
        segments,
        gc_tag,
    })
}

/// Build a fragment that exposes k-mer-safe reference segments.
///
/// The function pairs two reads, honours the requested indel handling strategy, optionally merges
/// the inter-mate gap, and trims user-defined offsets from both fragment ends. Only bases that are
/// simultaneously reliable for k-mer enumeration survive into the final segment list.
///
/// Parameters:
/// - `a`: First read to pair.
/// - `b`: Second read to pair.
/// - `indel_mode`: Controls whether reads containing insertions/deletions are skipped or flattened.
/// - `include_inter_mate_gap`: When true, the span between the mates is folded into the neighbouring
///   segment(s) whenever both sides provide sequence.
/// - `end_offset`: Number of bases to trim from each fragment end.
///
/// Returns `None` when the read pair fails orientation checks, is filtered out by `indel_mode`, or
/// trimming removes all reference sequence.
pub(crate) fn collect_fragment_with_kmer_segments(
    a: &KmerSegmentedReadInfo,
    b: &KmerSegmentedReadInfo,
    indel_mode: IndelMode,
    include_inter_mate_gap: bool,
    end_offset: u32,
) -> Option<FragmentWithKmerSegments> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }
    let gc_tag = combine_gc_tag_values(&forward.gc_tag, &reverse.gc_tag);

    if matches!(indel_mode, IndelMode::Skip) && (forward.has_indel() || reverse.has_indel()) {
        return None;
    }

    let span_start = forward.start();
    let span_end = reverse.end();

    let treat_as_flat =
        matches!(indel_mode, IndelMode::Ignore) || (!forward.has_indel() && !reverse.has_indel());

    if treat_as_flat {
        return collect_flat_fragment(
            forward,
            reverse,
            include_inter_mate_gap,
            end_offset,
            span_start,
            span_end,
            gc_tag,
        );
    }

    // Build absolute segments for each mate. Indel adjustments keep the fine-grained spans so we can
    // exclude bases touched by both reads' indels when requested.
    let mut forward_segments = forward.absolute_segments();
    let mut reverse_segments = reverse.absolute_segments();

    // Optionally bridge the mate gap. We extend into the gap only when the neighbouring edge is
    // backed by reference sequence (no trailing/leading insertion). If both sides carry indels, the
    // gap becomes its own segment so callers can still opt-in to counting within it.
    if include_inter_mate_gap && forward.end() < reverse.start() {
        let gap_interval = Interval::new(forward.end(), reverse.start()).ok()?;
        let mut left_extended = false;
        if let Some(last) = forward_segments.last_mut()
            && last.end() == forward.end()
            && !forward.trailing_insertion
        {
            *last = Interval::new(last.start(), gap_interval.end()).ok()?;
            left_extended = true;
        }

        let mut right_extended = false;
        if let Some(first) = reverse_segments.first_mut()
            && first.start() == reverse.start()
            && !reverse.leading_insertion
        {
            *first = Interval::new(gap_interval.start(), first.end()).ok()?;
            right_extended = true;
        }

        if !left_extended && !right_extended {
            forward_segments.push(gap_interval);
        }
    }

    if forward_segments.is_empty() && reverse_segments.is_empty() {
        return None;
    }

    forward_segments.sort_unstable_by_key(|segment| segment.start());
    reverse_segments.sort_unstable_by_key(|segment| segment.start());

    // Trim fixed offsets from both ends so k-mer contexts avoid edge artifacts.
    // `end_offset` is expected to be small (defaults to 0), so most spans dwarf it.
    // When a fragment is shorter than 2 * end_offset, trimming would invert the span.
    // We thus guard by collapsing `trim_end` to `span_start` in that case so the
    // subsequent check returns None.
    let trim_start = span_start.saturating_add(end_offset);
    let trim_end = if span_end > end_offset {
        span_end - end_offset
    } else {
        span_start
    };
    let trim_window = Interval::new(trim_start, trim_end).ok()?;

    let mut candidates: Vec<Interval<u32>> = forward_segments
        .into_iter()
        .chain(reverse_segments)
        .filter_map(|segment| segment.clip_to(trim_window))
        .collect();

    if candidates.is_empty() {
        return None;
    }

    candidates.sort_unstable_by_key(|segment| segment.start());

    // Keep bases covered by either mate. Touching spans stay separate so we
    // keep hard boundaries around insertions.
    let merged = merge_sorted_intervals(candidates, TouchingMergePolicy::KeepTouchingSeparate);

    if merged.is_empty() {
        return None;
    }

    let mut segments: SmallVec<[Interval<u32>; 12]> = SmallVec::with_capacity(merged.len());
    segments.extend(merged);

    if segments.is_empty() {
        return None;
    }

    Some(FragmentWithKmerSegments {
        interval: Interval::new(span_start, span_end).ok()?,
        segments,
        gc_tag,
    })
}

/// Build a fragment with k-mer-safe segments from a single read (unpaired input).
pub(crate) fn collect_fragment_with_kmer_segments_from_single_read(
    read: &KmerSegmentedReadInfo,
    indel_mode: IndelMode,
    end_offset: u32,
) -> Option<FragmentWithKmerSegments> {
    if matches!(indel_mode, IndelMode::Skip) && read.has_indel() {
        return None;
    }

    let span_start = read.start();
    let span_end = read.end();

    // Trim fixed offsets from both ends so k-mer contexts avoid edge artifacts.
    // `end_offset` is expected to be small (defaults to 0), so most spans dwarf it.
    // When a fragment is shorter than 2 * end_offset, trimming would invert the span.
    // We thus guard by collapsing `trim_end` to `span_start` in that case so the
    // subsequent check returns None.
    let trim_start = span_start.saturating_add(end_offset);
    let trim_end = if span_end > end_offset {
        span_end - end_offset
    } else {
        span_start
    };
    let trim_window = Interval::new(trim_start, trim_end).ok()?;

    let treat_as_flat = matches!(indel_mode, IndelMode::Ignore) || !read.has_indel();

    let mut segments: SmallVec<[Interval<u32>; 12]> = SmallVec::new();

    if treat_as_flat {
        segments.push(trim_window);
    } else {
        let mut candidates: Vec<Interval<u32>> = read
            .absolute_segments()
            .into_iter()
            .filter_map(|segment| segment.clip_to(trim_window))
            .collect();

        if candidates.is_empty() {
            return None;
        }

        candidates.sort_unstable_by_key(|segment| segment.start());

        let merged = merge_sorted_intervals(candidates, TouchingMergePolicy::KeepTouchingSeparate);

        if merged.is_empty() {
            return None;
        }

        segments.extend(merged);
    }

    if segments.is_empty() {
        return None;
    }

    Some(FragmentWithKmerSegments {
        interval: Interval::new(span_start, span_end).ok()?,
        segments,
        gc_tag: read.gc_tag,
    })
}

#[cfg(test)]
mod tests {
    include!("segment_kmer_fragment_tests.rs");
}
