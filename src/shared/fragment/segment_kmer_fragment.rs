use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, Record};
use smallvec::SmallVec;

use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::gc_tag::{GcTagValue, combine_gc_tag_values, read_gc_tag_from_record};
use crate::shared::indel_mode::IndelMode;

/// Represents a fragment together with the reference segments that are safe for k-mer analysis.
#[derive(Debug, Clone)]
pub struct FragmentWithKmerSegments {
    pub tid: i32,
    pub start: u32,
    pub end: u32,
    pub segments: SmallVec<[(u32, u32); 12]>,
    pub gc_tag: GcTagValue,
}

impl FragmentWithKmerSegments {
    #[inline]
    pub fn len(&self) -> u32 {
        self.end.saturating_sub(self.start)
    }

    #[inline]
    pub fn total_segment_len(&self) -> u32 {
        self.segments
            .iter()
            .fold(0u32, |acc, (s, e)| acc + e.saturating_sub(*s))
    }
}

/// Captures per-read segment data together with indel information that influences k-mer counting.
#[derive(Debug, Clone)]
pub struct KmerSegmentedReadInfo {
    pub tid: i32,
    pub pos: u32,
    pub end: u32,
    pub is_reverse: bool,
    pub has_insertion: bool,
    pub has_deletion: bool,
    pub leading_insertion: bool,
    pub trailing_insertion: bool,
    pub ref_mapped_segments: Vec<(u32, u32)>,
    pub gc_tag: GcTagValue,
}

impl KmerSegmentedReadInfo {
    #[inline]
    pub fn has_indel(&self) -> bool {
        self.has_insertion || self.has_deletion
    }

    #[inline]
    pub fn absolute_segments(&self) -> Vec<(u32, u32)> {
        if self.ref_mapped_segments.is_empty() {
            vec![(self.pos, self.end)]
        } else {
            self.ref_mapped_segments
                .iter()
                .map(|(off, len)| {
                    let start = self.pos.saturating_add(*off);
                    let end = start.saturating_add(*len);
                    (start, end)
                })
                .collect()
        }
    }
}

impl KmerSegmentedReadInfo {
    /// Build read metadata, optionally collecting reference segments for indel-aware counting.
    pub fn from_record(r: &Record, capture_segments: bool, gc_tag: Option<&[u8]>) -> Self {
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

        KmerSegmentedReadInfo {
            tid: r.tid(),
            pos: r.pos() as u32,
            end: r.reference_end() as u32,
            is_reverse: r.is_reverse(),
            has_insertion,
            has_deletion,
            leading_insertion,
            trailing_insertion,
            ref_mapped_segments,
            gc_tag: gc_tag_value,
        }
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
        self.pos
    }
}

/// Clip an interval to the `[trim_start, trim_end)` window.
#[inline]
fn clip_interval(start: u32, end: u32, trim_start: u32, trim_end: u32) -> Option<(u32, u32)> {
    let clipped_start = start.max(trim_start);
    let clipped_end = end.min(trim_end);
    if clipped_end > clipped_start {
        Some((clipped_start, clipped_end))
    } else {
        None
    }
}

/// Append a segment, coalescing it with the previous one when they overlap or touch.
#[inline]
fn push_merged(segments: &mut Vec<(u32, u32)>, segment: (u32, u32)) {
    if segment.1 <= segment.0 {
        return;
    }

    if let Some(last) = segments.last_mut() {
        if last.1 > segment.0 {
            if segment.1 > last.1 {
                last.1 = segment.1;
            }
            return;
        }
    }

    // Segments that just touch get appended as-is so we keep hard boundaries around insertions.
    segments.push(segment);
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
    gc_tag: GcTagValue,
) -> Option<FragmentWithKmerSegments> {
    let trim_start = span_start.saturating_add(end_offset);
    let trim_end = if span_end > end_offset {
        span_end - end_offset
    } else {
        span_start
    };
    if trim_start >= trim_end {
        return None;
    }

    let forward_segment = clip_interval(forward.pos, forward.end, trim_start, trim_end);
    let reverse_segment = clip_interval(reverse.pos, reverse.end, trim_start, trim_end);

    let mates_overlap_or_touch = forward.end >= reverse.pos;
    let mut segments: SmallVec<[(u32, u32); 12]> = SmallVec::new();

    if include_inter_mate_gap || mates_overlap_or_touch {
        segments.push((trim_start, trim_end));
    } else {
        if let Some((fs, fe)) = forward_segment {
            segments.push((fs, fe));
        }

        if let Some((rs, re)) = reverse_segment {
            if let Some(last) = segments.last_mut() {
                if last.1 >= rs {
                    if re > last.1 {
                        last.1 = re;
                    }
                } else {
                    segments.push((rs, re));
                }
            } else {
                segments.push((rs, re));
            }
        }
    }

    if segments.is_empty() {
        return None;
    }

    Some(FragmentWithKmerSegments {
        tid: forward.tid,
        start: span_start,
        end: span_end,
        segments,
        gc_tag,
    })
}

/// Build a fragment that exposes k-mer safe reference segments.
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
pub fn collect_fragment_with_kmer_segments(
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

    let span_start = forward.pos;
    let span_end = reverse.end;
    if span_start >= span_end {
        return None;
    }

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
    if include_inter_mate_gap && forward.end < reverse.pos {
        let mut left_extended = false;
        if let Some(last) = forward_segments.last_mut() {
            if last.1 == forward.end && !forward.trailing_insertion {
                last.1 = reverse.pos;
                left_extended = true;
            }
        }

        let mut right_extended = false;
        if let Some(first) = reverse_segments.first_mut() {
            if first.0 == reverse.pos && !reverse.leading_insertion {
                first.0 = forward.end;
                right_extended = true;
            }
        }

        if !left_extended && !right_extended {
            forward_segments.push((forward.end, reverse.pos));
        }
    }

    if forward_segments.is_empty() && reverse_segments.is_empty() {
        return None;
    }

    forward_segments.sort_unstable_by_key(|&(s, _)| s);
    reverse_segments.sort_unstable_by_key(|&(s, _)| s);

    let trim_start = span_start.saturating_add(end_offset);
    let trim_end = if span_end > end_offset {
        span_end - end_offset
    } else {
        span_start
    };
    if trim_start >= trim_end {
        return None;
    }

    let mut candidates: Vec<(u32, u32)> = forward_segments
        .into_iter()
        .chain(reverse_segments.into_iter())
        .filter_map(|(s, e)| clip_interval(s, e, trim_start, trim_end))
        .collect();

    if candidates.is_empty() {
        return None;
    }

    candidates.sort_unstable_by_key(|&(s, _)| s);

    let mut merged: Vec<(u32, u32)> = Vec::new();
    for seg in candidates.into_iter() {
        // Keep bases covered by either mate; overlapping spans collapse via `push_merged` so we only
        // exclude reference positions when both reads align an indel there.
        push_merged(&mut merged, seg);
    }

    if merged.is_empty() {
        return None;
    }

    let mut segments: SmallVec<[(u32, u32); 12]> = SmallVec::with_capacity(merged.len());
    segments.extend(merged.into_iter());

    if segments.is_empty() {
        return None;
    }

    Some(FragmentWithKmerSegments {
        tid: forward.tid,
        start: span_start,
        end: span_end,
        segments,
        gc_tag,
    })
}
