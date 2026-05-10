use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, Record};
use smallvec::SmallVec;

use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    Fragment, PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::gc_tag::{GCTagValue, combine_gc_tag_values, read_gc_tag_from_record};
use crate::shared::interval::{Interval, TouchingMergePolicy, merge_sorted_intervals};

/// Fragment that may carry explicit reference-coverage segments
///
/// If `segments` is `None`, use the plain fragment interval `[start, end)`
/// If `segments` is `Some`, use those `[start, end)` segments instead
#[derive(Debug, Clone)]
pub struct FragmentWithSegments {
    pub tid: i32,
    pub interval: Interval<u32>, // forward.start .. reverse.end
    pub segments: Option<SmallVec<[Interval<u32>; 12]>>,
    pub gc_tag: GCTagValue,
}

impl FragmentWithSegments {
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Length of the fragment (end - start).
    #[inline]
    pub fn len(&self) -> u32 {
        self.interval.len()
    }
}

impl From<Fragment> for FragmentWithSegments {
    fn from(f: Fragment) -> Self {
        FragmentWithSegments {
            tid: f.tid,
            interval: f.interval,
            segments: None,
            gc_tag: GCTagValue::default(),
        }
    }
}

/// Compact per-read metadata plus optional mapped-reference segments
///
/// Stores only what we need to assemble a fragment without keeping whole BAM records.
/// If the read's CIGAR contains reference gaps (`D` or `N`), we also store the
/// read's **mapped reference segments** as relative pairs `[offset_from_pos, len]`
/// for ref+query consuming ops (`M`, `=`, `X`). Otherwise `ref_mapped_segments` is empty.
///
/// Notes
/// -----
/// - `interval` is the read's aligned reference span
/// - `ref_mapped_segments` elements are relative to `interval.start()`
/// - Adjacent segments separated only by non-reference ops are merged
#[derive(Debug, Clone)]
pub struct SegmentedReadInfo {
    pub tid: i32,                // Contig id
    pub interval: Interval<u32>, // Aligned reference span [start: pos(), end: reference_end())
    pub is_reverse: bool,
    pub has_ref_gap: bool,                    // True if any D/N present
    pub max_ref_gap: u32,                     // Longest single D/N length (0 if none)
    pub ref_mapped_segments: Vec<(u32, u32)>, // Relative segments: (offset_from_pos, len)
    pub gc_tag: GCTagValue,
}

impl SegmentedReadInfo {
    /// Return the read's inclusive start on the reference.
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Return the read's exclusive end on the reference.
    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Return the read's aligned reference span `[pos, end)`.
    #[inline]
    pub fn aligned_interval(&self) -> Interval<u32> {
        self.interval
    }

    #[inline]
    pub fn from_record_with_gc_tag(r: &Record, gc_tag: Option<&[u8]>) -> Result<Self> {
        // Detect any D/N and track max gap length
        let mut has_ref_gap = false;
        let mut max_gap: u32 = 0;
        for c in r.cigar().iter() {
            match *c {
                Cigar::Del(l) | Cigar::RefSkip(l) => {
                    has_ref_gap = true;
                    if l > max_gap {
                        max_gap = l;
                    }
                }
                _ => {}
            }
        }

        // Build relative ref+query segments only if needed
        let mut ref_mapped_segments: Vec<(u32, u32)> = Vec::new();
        if has_ref_gap {
            let mut ref_off: u32 = 0;
            let mut seg_start: Option<u32> = None;

            for c in r.cigar().iter() {
                match *c {
                    // Consume ref+query -> extend or start a segment
                    Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                        if seg_start.is_none() {
                            seg_start = Some(ref_off);
                        }
                        ref_off += l;
                    }
                    // Consume ref only (gap) -> close open segment and advance
                    Cigar::Del(l) | Cigar::RefSkip(l) => {
                        if let Some(s) = seg_start.take() {
                            let len = ref_off.saturating_sub(s);
                            if len > 0 {
                                ref_mapped_segments.push((s, len));
                            }
                        }
                        ref_off += l;
                    }
                    // Do not consume ref -> no advance
                    Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
            }
            // Close trailing segment if open
            if let Some(s) = seg_start.take() {
                let len = ref_off.saturating_sub(s);
                if len > 0 {
                    ref_mapped_segments.push((s, len));
                }
            }
            // Merge overlaps/adjacency on reference
            if ref_mapped_segments.len() > 1 {
                ref_mapped_segments.sort_unstable_by_key(|&(off, _)| off);
                let mut merged: Vec<(u32, u32)> = Vec::with_capacity(ref_mapped_segments.len());
                for (off, len) in ref_mapped_segments.drain(..) {
                    if let Some(last) = merged.last_mut() {
                        let last_end = last.0 + last.1;
                        if off <= last_end {
                            let end = (off + len).max(last_end);
                            last.1 = end - last.0;
                            continue;
                        }
                    }
                    merged.push((off, len));
                }
                ref_mapped_segments = merged;
            }
        }

        let gc_tag_value = gc_tag
            .map(|tag| read_gc_tag_from_record(r, tag))
            .unwrap_or_default();

        Ok(SegmentedReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            has_ref_gap,
            max_ref_gap: max_gap,
            ref_mapped_segments,
            gc_tag: gc_tag_value,
        })
    }
}

impl PairOrientable for SegmentedReadInfo {
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

/// Build a `FragmentWithSegments` from two `SegmentedReadInfo` instances
///
/// Summary
/// -------
/// Returns a fragment spanning `[forward.pos, reverse.reference_end)` and, when either read
/// has a sufficiently large reference gap, attaches explicit mapped-reference
/// segments so downstream coverage respects true deletions.
///
/// Parameters
/// ----------
/// - a: First read
/// - b: Mate read
/// - trigger_min_gap_bp: Minimum D/N length in either read to trigger segment mode
/// - include_inter_mate_gap: Count the [forward.end, reverse.pos) gap
///   (when reads don't overlap) as part of the fragment
///
/// Returns
/// -------
/// - frag: FragmentWithSegments covering `[forward.pos, reverse.reference_end)`
///   With `segments = None` when no triggering gap is present
pub fn collect_fragment_with_segments(
    a: &SegmentedReadInfo,
    b: &SegmentedReadInfo,
    trigger_min_gap_bp: u32,
    include_inter_mate_gap: bool,
) -> Option<FragmentWithSegments> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }

    let fragment_interval = Interval::new(forward.start(), reverse.end()).ok()?;
    let gc_tag = combine_gc_tag_values(&forward.gc_tag, &reverse.gc_tag);

    // Decide if we switch to segments
    let trigger = (forward.has_ref_gap && forward.max_ref_gap >= trigger_min_gap_bp)
        || (reverse.has_ref_gap && reverse.max_ref_gap >= trigger_min_gap_bp);

    // If no trigger and user wants full fragment counting, return the plain fragment interval
    if !trigger && include_inter_mate_gap {
        return Some(FragmentWithSegments {
            tid: forward.tid,
            interval: fragment_interval,
            segments: None,
            gc_tag,
        });
    }

    // Build absolute segments to honor either:
    // - Triggered ref gaps, and optionally add the inter-mate gap
    // - Or, when not triggered and include_inter_mate_gap == false (the +2),
    //   exclude the inter-mate gap by using only per-read spans
    let mut abs: Vec<Interval<u32>> = Vec::with_capacity(
        2 + forward
            .ref_mapped_segments
            .len()
            .saturating_add(reverse.ref_mapped_segments.len()),
    );

    // Expand forward read's relative ref-mapped segments to absolute coordinates
    //
    // Each stored tuple is (offset_from_pos, len) measured on the reference
    // Add `pos` to get absolute [start, end) on the chromosome
    // If the list is empty (no gaps worth storing), fall back to the read's aligned interval [pos, end)
    if !forward.ref_mapped_segments.is_empty() {
        for (off, len) in &forward.ref_mapped_segments {
            let segment_start = forward.start().saturating_add(*off);
            let segment_end = segment_start.saturating_add(*len);
            let Ok(segment) = Interval::new(segment_start, segment_end) else {
                continue;
            };
            abs.push(segment);
        }
    } else {
        abs.push(forward.aligned_interval());
    }

    // Same expansion for the reverse read
    if !reverse.ref_mapped_segments.is_empty() {
        for (off, len) in &reverse.ref_mapped_segments {
            let segment_start = reverse.start().saturating_add(*off);
            let segment_end = segment_start.saturating_add(*len);
            let Ok(segment) = Interval::new(segment_start, segment_end) else {
                continue;
            };
            abs.push(segment);
        }
    } else {
        abs.push(reverse.aligned_interval());
    }

    // Optionally include the fragment insert between mates
    //
    // Rationale: For fragment coverage, some users want the full molecule counted
    // across the unsequenced insert between the reads. That gap is the reference
    // span from the end of the forward read to the start of the reverse read
    // (when they do not overlap). If reads overlap, there is no gap to add
    //
    // Note: When !trigger and include_inter_mate_gap == false we intentionally do not add the gap
    if trigger && include_inter_mate_gap && forward.end() < reverse.start() {
        abs.push(Interval::new(forward.end(), reverse.start()).ok()?);
    }

    if abs.is_empty() {
        // Fallback to the plain fragment interval
        return Some(FragmentWithSegments {
            tid: forward.tid,
            interval: fragment_interval,
            segments: None,
            gc_tag,
        });
    }

    abs.sort_unstable_by_key(|segment| segment.start());

    // Merge and clip to the fragment interval
    let merged = merge_sorted_intervals(
        abs.into_iter()
            .filter_map(|segment| segment.clip_to(fragment_interval))
            .collect(),
        TouchingMergePolicy::MergeTouching,
    );

    let segments = if merged.is_empty() {
        None
    } else {
        Some(merged.into())
    };

    Some(FragmentWithSegments {
        tid: forward.tid,
        interval: fragment_interval,
        segments,
        gc_tag,
    })
}

/// Build a fragment from a single segmented read (unpaired input).
pub fn collect_fragment_with_segments_from_single_read(
    read: &SegmentedReadInfo,
    trigger_min_gap_bp: u32,
) -> Option<FragmentWithSegments> {
    let fragment_interval = read.aligned_interval();

    // Decide if we switch to segments based on reference gaps
    let trigger = read.has_ref_gap && read.max_ref_gap >= trigger_min_gap_bp;

    // If no trigger, return the plain fragment interval
    if !trigger {
        return Some(FragmentWithSegments {
            tid: read.tid,
            interval: fragment_interval,
            segments: None,
            gc_tag: read.gc_tag,
        });
    }

    // Expand reference-mapped segments to absolute coordinates
    //
    // Each stored tuple is (offset_from_pos, len) measured on the reference
    // Add `pos` to get absolute [start, end) on the chromosome
    let mut abs: Vec<Interval<u32>> = Vec::with_capacity(read.ref_mapped_segments.len());
    if !read.ref_mapped_segments.is_empty() {
        for (off, len) in &read.ref_mapped_segments {
            let segment_start = read.start().saturating_add(*off);
            let segment_end = segment_start.saturating_add(*len);
            let Ok(segment) = Interval::new(segment_start, segment_end) else {
                continue;
            };
            abs.push(segment);
        }
    }

    if abs.is_empty() {
        return Some(FragmentWithSegments {
            tid: read.tid,
            interval: fragment_interval,
            segments: None,
            gc_tag: read.gc_tag,
        });
    }

    // Segments are already merged and sorted in `SegmentedReadInfo::from_record_with_gc_tag`
    // so we can attach them directly. Keep a light validity check only.
    let segments = if abs.is_empty() {
        None
    } else {
        let mut v = SmallVec::with_capacity(abs.len());
        for segment in abs {
            if let Some(segment) = segment.clip_to(fragment_interval) {
                v.push(segment);
            }
        }
        if v.is_empty() { None } else { Some(v) }
    };

    Some(FragmentWithSegments {
        tid: read.tid,
        interval: fragment_interval,
        segments,
        gc_tag: read.gc_tag,
    })
}
