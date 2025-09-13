use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::{Cigar, Record};
use smallvec::SmallVec;

use crate::utils::fragment::minimal_fragment::{
    Fragment, PairOrientable, oriented_pair_from_read_info,
};

/// Fragment that may carry explicit reference-coverage segments
///
/// If `segments` is `None`, use the plain fragment span `[start, end)`
/// If `segments` is `Some`, use those `[start, end)` segments instead
#[derive(Debug, Clone)]
pub struct FragmentWithSegments {
    pub tid: i32,
    pub start: u32, // forward.start
    pub end: u32,   // reverse.end (end-exclusive)
    pub segments: Option<SmallVec<[(u32, u32); 12]>>,
}

impl From<Fragment> for FragmentWithSegments {
    fn from(f: Fragment) -> Self {
        FragmentWithSegments {
            tid: f.tid,
            start: f.start,
            end: f.end,
            segments: None,
        }
    }
}

/// Compact per-read metadata plus optional mapped-reference segments
///
/// Stores only what we need to assemble a fragment without keeping whole BAM records.
/// If the read’s CIGAR contains reference gaps (`D` or `N`), we also store the
/// read’s **mapped reference segments** as relative pairs `[offset_from_pos, len]`
/// for ref+query consuming ops (`M`, `=`, `X`). Otherwise `ref_mapped_segments` is empty.
///
/// Notes
/// -----
/// - `pos` and `end` are the read’s aligned reference coordinates
/// - `ref_mapped_segments` elements are relative to `pos`
/// - Adjacent segments separated only by non-reference ops are merged
#[derive(Debug, Clone)]
pub struct SegmentedReadInfo {
    pub tid: i32, // Contig id
    pub pos: u32, // 0-based leftmost aligned ref pos
    pub end: u32, // 0-based exclusive rightmost aligned ref pos
    pub is_reverse: bool,
    pub has_ref_gap: bool,                    // True if any D/N present
    pub max_ref_gap: u32,                     // Longest single D/N length (0 if none)
    pub ref_mapped_segments: Vec<(u32, u32)>, // Relative segments: (offset_from_pos, len)
}

impl From<&Record> for SegmentedReadInfo {
    #[inline]
    fn from(r: &Record) -> Self {
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

        SegmentedReadInfo {
            tid: r.tid(),
            pos: r.pos() as u32,
            end: r.reference_end() as u32,
            is_reverse: r.is_reverse(),
            has_ref_gap,
            max_ref_gap: max_gap,
            ref_mapped_segments,
        }
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
}

/// Build a `FragmentWithSegments` from two `SegmentedReadInfo` instances
///
/// Summary
/// -------
/// Returns a fragment spanning `[forward.pos, reverse.end)` and, when either read
/// has a sufficiently large reference gap, attaches explicit mapped-reference
/// segments so downstream coverage respects true deletions.
///
/// Parameters
/// ----------
/// - a: First read
/// - b: Mate read
/// - trigger_min_gap_bp: Minimum D/N length in either read to trigger segment mode
/// - include_inter_mate_gap: Count the [fwd.end, rev.pos) gap
///   (when reads don't overlap) as part of the fragment
///
/// Returns
/// -------
/// - frag: FragmentWithSegments covering `[forward.pos, reverse.end)`
///   With `segments = None` when no triggering gap is present
pub fn collect_fragment_with_segments(
    a: &SegmentedReadInfo,
    b: &SegmentedReadInfo,
    trigger_min_gap_bp: u32,
    include_inter_mate_gap: bool,
) -> Option<FragmentWithSegments> {
    let (fwd, rev) = oriented_pair_from_read_info(a, b)?;
    if rev.end <= fwd.pos {
        return None;
    }

    let span_start = fwd.pos;
    let span_end = rev.end;

    // Decide if we switch to segments
    let trigger = (fwd.has_ref_gap && fwd.max_ref_gap >= trigger_min_gap_bp)
        || (rev.has_ref_gap && rev.max_ref_gap >= trigger_min_gap_bp);

    // If no trigger and user wants full fragment counting, return the plain span
    if !trigger && include_inter_mate_gap {
        return Some(FragmentWithSegments {
            tid: fwd.tid,
            start: span_start,
            end: span_end,
            segments: None,
        });
    }

    // Build absolute segments to honor either:
    // - Triggered ref gaps, and optionally add the inter-mate gap
    // - Or, when not triggered and include_inter_mate_gap == false (the +2),
    //   exclude the inter-mate gap by using only per-read spans
    let mut abs: Vec<(u32, u32)> = Vec::with_capacity(
        2 + fwd
            .ref_mapped_segments
            .len()
            .saturating_add(rev.ref_mapped_segments.len()),
    );

    // Expand forward read's relative ref-mapped segments to absolute genome coords
    //
    // Each stored tuple is (offset_from_pos, len) measured on the reference
    // Add `pos` to get absolute [start, end) on the chromosome
    // If the list is empty (no gaps worth storing), fall back to the read's aligned span [pos, end)
    if !fwd.ref_mapped_segments.is_empty() {
        for (off, len) in &fwd.ref_mapped_segments {
            let s = fwd.pos.saturating_add(*off);
            let e = s.saturating_add(*len);
            abs.push((s, e));
        }
    } else {
        abs.push((fwd.pos, fwd.end));
    }

    // Same expansion for the reverse read
    if !rev.ref_mapped_segments.is_empty() {
        for (off, len) in &rev.ref_mapped_segments {
            let s = rev.pos.saturating_add(*off);
            let e = s.saturating_add(*len);
            abs.push((s, e));
        }
    } else {
        abs.push((rev.pos, rev.end));
    }

    // Optionally include the fragment insert between mates
    //
    // Rationale: For fragment coverage, some users want the full molecule counted
    // across the unsequenced insert between the reads. That gap is the reference
    // span from the end of the forward read to the start of the reverse read
    // (when they do not overlap). If reads overlap, there is no gap to add
    //
    // Note: When !trigger and include_inter_mate_gap == false we intentionally do not add the gap
    if trigger && include_inter_mate_gap {
        if fwd.end < rev.pos {
            abs.push((fwd.end, rev.pos));
        }
    }

    if abs.is_empty() {
        // Fallback to plain span
        return Some(FragmentWithSegments {
            tid: fwd.tid,
            start: span_start,
            end: span_end,
            segments: None,
        });
    }

    abs.sort_unstable_by_key(|&(s, _)| s);

    // Merge and clip to fragment span
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(abs.len());
    for (mut s, mut e) in abs {
        // Check validity of segment
        if s >= e {
            continue;
        }
        if e <= span_start || s >= span_end {
            continue;
        }

        // Clip to span
        if s < span_start {
            s = span_start;
        }
        if e > span_end {
            e = span_end;
        }

        // Merge overlapping segments by increasing the previous segment
        if let Some(last) = merged.last_mut() {
            if s <= last.1 {
                if e > last.1 {
                    last.1 = e;
                }
                continue;
            }
        }
        merged.push((s, e));
    }

    let segments = if merged.is_empty() {
        None
    } else {
        Some(merged.into())
    };

    Some(FragmentWithSegments {
        tid: fwd.tid,
        start: span_start,
        end: span_end,
        segments,
    })
}
