use crate::{
    Result,
    commands::ends::config_structs::{ClipStrategy, KmerSource},
    shared::{
        fragment::minimal_fragment::{
            PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
        },
        gc_tag::{GcTagValue, combine_gc_tag_values, read_gc_tag_from_record},
        indel_mode::IndelMotifFilterPolicy,
        interval::Interval,
    },
};
use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Cigar, Record},
};

/// One counted fragment end after clip handling has been resolved.
#[derive(Debug, Clone)]
pub struct ResolvedFragmentEnd {
    /// Fragment boundary used by the selected clip strategy.
    pub boundary_pos: u32,
    /// First `k_inside` bases adjacent to this end in BAM/reference storage orientation.
    pub inside_bases: Vec<u8>,
}

/// Fragment payload for the `ends` command.
#[derive(Debug, Clone)]
pub struct FragmentWithEnds {
    pub tid: i32,
    /// Aligned fragment interval used for length and GC-related fragment geometry.
    pub interval: Interval<u32>,
    /// Boundary-adjusted interval used for window assignment when clip strategy changes the ends.
    pub assignment_interval: Interval<u32>,
    pub gc_tag: GcTagValue,
    pub left_end: Option<ResolvedFragmentEnd>,
    pub right_end: Option<ResolvedFragmentEnd>,
}

impl FragmentWithEnds {
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }

    #[inline]
    pub fn len(&self) -> u32 {
        self.interval.len()
    }

    #[inline]
    pub fn assignment_len(&self) -> u32 {
        self.assignment_interval.len()
    }

    #[inline]
    pub fn assignment_start(&self) -> u32 {
        self.assignment_interval.start()
    }

    #[inline]
    pub fn assignment_end(&self) -> u32 {
        self.assignment_interval.end()
    }
}

/// Compact per-read data needed to assemble `FragmentWithEnds`.
#[derive(Debug, Clone)]
pub struct EndReadInfo {
    pub tid: i32,
    pub interval: Interval<u32>,
    pub is_reverse: bool,
    pub left_soft_clip_bp: u32,
    pub right_soft_clip_bp: u32,
    pub left_motif_has_indels: bool,
    pub right_motif_has_indels: bool,
    pub has_hard_clip: bool,
    pub seq: Vec<u8>,
    pub gc_tag: GcTagValue,
}

impl EndReadInfo {
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }

    #[inline]
    pub fn aligned_interval(&self) -> Interval<u32> {
        self.interval
    }

    pub fn from_record_with_gc_tag(
        r: &Record,
        gc_tag: Option<&[u8]>,
        clip_strategy: ClipStrategy,
        k_inside: usize,
    ) -> Result<Self> {
        let (left_soft_clip_bp, right_soft_clip_bp, has_hard_clip) = inspect_cigar_edges(r);
        let gc_tag_value = gc_tag
            .map(|tag| read_gc_tag_from_record(r, tag))
            .unwrap_or_default();

        Ok(EndReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            left_soft_clip_bp,
            right_soft_clip_bp,
            left_motif_has_indels: motif_has_indels(
                r,
                FragmentEndSide::Left,
                clip_strategy,
                k_inside,
                left_soft_clip_bp,
            ),
            right_motif_has_indels: motif_has_indels(
                r,
                FragmentEndSide::Right,
                clip_strategy,
                k_inside,
                right_soft_clip_bp,
            ),
            has_hard_clip,
            seq: r.seq().as_bytes(),
            gc_tag: gc_tag_value,
        })
    }
}

impl PairOrientable for EndReadInfo {
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

#[derive(Debug, Clone, Copy)]
enum FragmentEndSide {
    Left,
    Right,
}

#[derive(Debug, Clone)]
enum ResolvedEndOutcome {
    /// Abort the whole fragment, e.g. `skip-affected-fragment`.
    DropFragment,
    /// Skip this end's motif but still use its clip-adjusted assignment boundary.
    SkipEndKeepAssignmentBoundary { assignment_boundary_pos: u32 },
    /// Skip this end's motif and fall back to the aligned boundary for assignment.
    SkipEndDropAssignmentBoundary,
    /// Keep both fragment geometry and this end's motif.
    KeepEnd {
        assignment_boundary_pos: u32,
        end: ResolvedFragmentEnd,
    },
}

/// Build a `FragmentWithEnds` from two per-read summaries.
pub fn collect_fragment_with_ends(
    a: &EndReadInfo,
    b: &EndReadInfo,
    clip_strategy: ClipStrategy,
    source_inside: KmerSource,
    indel_filter: IndelMotifFilterPolicy,
    k_inside: usize,
    max_soft_clips: Option<u32>,
) -> Option<FragmentWithEnds> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }
    if forward.has_hard_clip || reverse.has_hard_clip {
        return None;
    }

    let (left_assignment_boundary_pos, left_end) = match resolve_fragment_end(
        forward,
        FragmentEndSide::Left,
        clip_strategy,
        source_inside,
        indel_filter,
        k_inside,
        max_soft_clips,
    ) {
        ResolvedEndOutcome::DropFragment => return None,
        ResolvedEndOutcome::SkipEndKeepAssignmentBoundary {
            assignment_boundary_pos,
        } => (assignment_boundary_pos, None),
        ResolvedEndOutcome::SkipEndDropAssignmentBoundary => (forward.start(), None),
        ResolvedEndOutcome::KeepEnd {
            assignment_boundary_pos,
            end,
        } => (assignment_boundary_pos, Some(end)),
    };
    let (right_assignment_boundary_pos, right_end) = match resolve_fragment_end(
        reverse,
        FragmentEndSide::Right,
        clip_strategy,
        source_inside,
        indel_filter,
        k_inside,
        max_soft_clips,
    ) {
        ResolvedEndOutcome::DropFragment => return None,
        ResolvedEndOutcome::SkipEndKeepAssignmentBoundary {
            assignment_boundary_pos,
        } => (assignment_boundary_pos, None),
        ResolvedEndOutcome::SkipEndDropAssignmentBoundary => (reverse.end(), None),
        ResolvedEndOutcome::KeepEnd {
            assignment_boundary_pos,
            end,
        } => (assignment_boundary_pos, Some(end)),
    };

    if left_end.is_none() && right_end.is_none() {
        return None;
    }

    let gc_tag = combine_gc_tag_values(&forward.gc_tag, &reverse.gc_tag);
    let interval = Interval::new(forward.start(), reverse.end()).ok()?;
    let assignment_interval =
        Interval::new(left_assignment_boundary_pos, right_assignment_boundary_pos).ok()?;

    Some(FragmentWithEnds {
        tid: forward.tid,
        interval,
        assignment_interval,
        gc_tag,
        left_end,
        right_end,
    })
}

/// Build a `FragmentWithEnds` from one read in `--reads-are-fragments` mode.
pub fn collect_fragment_with_ends_from_single_read(
    read: &EndReadInfo,
    clip_strategy: ClipStrategy,
    source_inside: KmerSource,
    indel_filter: IndelMotifFilterPolicy,
    k_inside: usize,
    max_soft_clips: Option<u32>,
) -> Option<FragmentWithEnds> {
    if read.has_hard_clip {
        return None;
    }

    let (left_assignment_boundary_pos, left_end) = match resolve_fragment_end(
        read,
        FragmentEndSide::Left,
        clip_strategy,
        source_inside,
        indel_filter,
        k_inside,
        max_soft_clips,
    ) {
        ResolvedEndOutcome::DropFragment => return None,
        ResolvedEndOutcome::SkipEndKeepAssignmentBoundary {
            assignment_boundary_pos,
        } => (assignment_boundary_pos, None),
        ResolvedEndOutcome::SkipEndDropAssignmentBoundary => (read.start(), None),
        ResolvedEndOutcome::KeepEnd {
            assignment_boundary_pos,
            end,
        } => (assignment_boundary_pos, Some(end)),
    };
    let (right_assignment_boundary_pos, right_end) = match resolve_fragment_end(
        read,
        FragmentEndSide::Right,
        clip_strategy,
        source_inside,
        indel_filter,
        k_inside,
        max_soft_clips,
    ) {
        ResolvedEndOutcome::DropFragment => return None,
        ResolvedEndOutcome::SkipEndKeepAssignmentBoundary {
            assignment_boundary_pos,
        } => (assignment_boundary_pos, None),
        ResolvedEndOutcome::SkipEndDropAssignmentBoundary => (read.end(), None),
        ResolvedEndOutcome::KeepEnd {
            assignment_boundary_pos,
            end,
        } => (assignment_boundary_pos, Some(end)),
    };

    if left_end.is_none() && right_end.is_none() {
        return None;
    }

    let interval = read.aligned_interval();
    let assignment_interval =
        Interval::new(left_assignment_boundary_pos, right_assignment_boundary_pos).ok()?;

    Some(FragmentWithEnds {
        tid: read.tid,
        interval,
        assignment_interval,
        gc_tag: read.gc_tag,
        left_end,
        right_end,
    })
}

/// Inspect terminal clipping in BAM storage orientation.
///
/// Valid SAM edge patterns are:
/// - `S...`
/// - `H...`
/// - `HS...`
/// - `...S`
/// - `...H`
/// - `...SH`
///
/// Soft clips contribute to the returned clip lengths.
/// Hard clips only contribute to the `has_hard_clip` flag.
fn inspect_cigar_edges(record: &Record) -> (u32, u32, bool) {
    let mut left_soft_clip_bp = 0;
    let mut right_soft_clip_bp = 0;
    let mut has_hard_clip = false;

    for op in record.cigar().iter() {
        match *op {
            Cigar::SoftClip(bp) => {
                left_soft_clip_bp += bp;
            }
            Cigar::HardClip(_) => {
                has_hard_clip = true;
            }
            _ => break,
        }
    }

    for op in record.cigar().iter().rev() {
        match *op {
            Cigar::SoftClip(bp) => {
                right_soft_clip_bp += bp;
            }
            Cigar::HardClip(_) => {
                has_hard_clip = true;
            }
            _ => break,
        }
    }

    (left_soft_clip_bp, right_soft_clip_bp, has_hard_clip)
}

/// Check whether the inside-fragment motif uses aligned bases across an indel.
fn motif_has_indels(
    record: &Record,
    end_side: FragmentEndSide,
    clip_strategy: ClipStrategy,
    k_inside: usize,
    soft_clip_bp: u32,
) -> bool {
    let aligned_bases_in_motif = match clip_strategy {
        ClipStrategy::Aligned | ClipStrategy::Drop => k_inside,
        ClipStrategy::Raw => k_inside.saturating_sub(soft_clip_bp as usize),
    };
    if aligned_bases_in_motif == 0 {
        return false;
    }

    let cigar_ops = record.cigar();
    let mut aligned_bases_seen = 0usize;

    match end_side {
        FragmentEndSide::Left => {
            for op in cigar_ops.iter() {
                if aligned_bases_seen >= aligned_bases_in_motif {
                    return false;
                }
                match *op {
                    Cigar::Ins(_) | Cigar::Del(_) | Cigar::RefSkip(_) => return true,
                    Cigar::Match(bp) | Cigar::Equal(bp) | Cigar::Diff(bp) => {
                        aligned_bases_seen = aligned_bases_seen.saturating_add(bp as usize);
                    }
                    Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
            }
        }
        FragmentEndSide::Right => {
            for op in cigar_ops.iter().rev() {
                if aligned_bases_seen >= aligned_bases_in_motif {
                    return false;
                }
                match *op {
                    Cigar::Ins(_) | Cigar::Del(_) | Cigar::RefSkip(_) => return true,
                    Cigar::Match(bp) | Cigar::Equal(bp) | Cigar::Diff(bp) => {
                        aligned_bases_seen = aligned_bases_seen.saturating_add(bp as usize);
                    }
                    Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
            }
        }
    }

    false
}

/// Resolve one fragment end according to the selected clip strategy.
///
/// Returns:
/// - `ResolvedEndOutcome::DropFragment` when the whole fragment must be discarded, e.g.
///   `skip-affected-fragment`
///
/// - `ResolvedEndOutcome::SkipEndKeepAssignmentBoundary` when this end's motif is skipped but the
///   end location should still use its clip-adjusted assignment boundary
///
/// - `ResolvedEndOutcome::SkipEndDropAssignmentBoundary` when this end is not trusted for motif
///   counting or assignment-boundary adjustment
///
/// - `ResolvedEndOutcome::KeepEnd` when both fragment geometry and this end's motif are kept
fn resolve_fragment_end(
    read: &EndReadInfo,
    end_side: FragmentEndSide,
    clip_strategy: ClipStrategy,
    source_inside: KmerSource,
    indel_filter: IndelMotifFilterPolicy,
    k_inside: usize,
    max_soft_clips: Option<u32>,
) -> ResolvedEndOutcome {
    let (aligned_boundary_pos, soft_clip_bp, motif_has_indels) = match end_side {
        FragmentEndSide::Left => (
            read.start(),
            read.left_soft_clip_bp,
            read.left_motif_has_indels,
        ),
        FragmentEndSide::Right => (
            read.end(),
            read.right_soft_clip_bp,
            read.right_motif_has_indels,
        ),
    };

    if max_soft_clips.is_some_and(|max_bp| soft_clip_bp > max_bp) {
        return ResolvedEndOutcome::SkipEndDropAssignmentBoundary;
    }

    // First decide whether indels invalidate just this end or the whole fragment.
    // Do not decide assignment-boundary behavior yet, because that depends on the clip strategy.
    if motif_has_indels && matches!(indel_filter, IndelMotifFilterPolicy::SkipAffectedFragment) {
        return ResolvedEndOutcome::DropFragment;
    }

    let skip_end_due_to_indels = motif_has_indels
        && match indel_filter {
            IndelMotifFilterPolicy::Auto => matches!(source_inside, KmerSource::Reference),
            IndelMotifFilterPolicy::SkipAffectedEnd
            | IndelMotifFilterPolicy::SkipAffectedFragment => true,
        };

    match clip_strategy {
        ClipStrategy::Aligned | ClipStrategy::Drop => {
            if matches!(clip_strategy, ClipStrategy::Drop) && soft_clip_bp > 0 {
                return ResolvedEndOutcome::SkipEndDropAssignmentBoundary;
            }
            if skip_end_due_to_indels {
                return ResolvedEndOutcome::SkipEndKeepAssignmentBoundary {
                    assignment_boundary_pos: aligned_boundary_pos,
                };
            }
            match build_resolved_end(
                read,
                end_side,
                clip_strategy,
                aligned_boundary_pos,
                k_inside,
            ) {
                Some(end) => ResolvedEndOutcome::KeepEnd {
                    assignment_boundary_pos: aligned_boundary_pos,
                    end,
                },
                None => ResolvedEndOutcome::SkipEndDropAssignmentBoundary,
            }
        }
        ClipStrategy::Raw => {
            let assignment_boundary_pos = match end_side {
                FragmentEndSide::Left => aligned_boundary_pos.saturating_sub(soft_clip_bp),
                FragmentEndSide::Right => aligned_boundary_pos.saturating_add(soft_clip_bp),
            };

            if skip_end_due_to_indels {
                return ResolvedEndOutcome::SkipEndKeepAssignmentBoundary {
                    assignment_boundary_pos,
                };
            }
            match build_resolved_end(
                read,
                end_side,
                clip_strategy,
                assignment_boundary_pos,
                k_inside,
            ) {
                Some(end) => ResolvedEndOutcome::KeepEnd {
                    assignment_boundary_pos,
                    end,
                },
                None => ResolvedEndOutcome::SkipEndDropAssignmentBoundary,
            }
        }
    }
}

fn build_resolved_end(
    read: &EndReadInfo,
    end_side: FragmentEndSide,
    clip_strategy: ClipStrategy,
    boundary_pos: u32,
    k_inside: usize,
) -> Option<ResolvedFragmentEnd> {
    let inside_bases = match end_side {
        FragmentEndSide::Left => extract_left_inside_bases(read, clip_strategy, k_inside)?,
        FragmentEndSide::Right => extract_right_inside_bases(read, clip_strategy, k_inside)?,
    };

    Some(ResolvedFragmentEnd {
        boundary_pos,
        inside_bases,
    })
}

fn extract_left_inside_bases(
    read: &EndReadInfo,
    clip_strategy: ClipStrategy,
    k_inside: usize,
) -> Option<Vec<u8>> {
    let start_idx = match clip_strategy {
        ClipStrategy::Aligned | ClipStrategy::Drop => read.left_soft_clip_bp as usize,
        ClipStrategy::Raw => 0,
    };
    let end_idx = start_idx.saturating_add(k_inside);
    if end_idx > read.seq.len() {
        return None;
    }
    Some(read.seq[start_idx..end_idx].to_vec())
}

fn extract_right_inside_bases(
    read: &EndReadInfo,
    clip_strategy: ClipStrategy,
    k_inside: usize,
) -> Option<Vec<u8>> {
    let end_idx = match clip_strategy {
        ClipStrategy::Aligned | ClipStrategy::Drop => read
            .seq
            .len()
            .checked_sub(read.right_soft_clip_bp as usize)?,
        ClipStrategy::Raw => read.seq.len(),
    };
    let start_idx = end_idx.checked_sub(k_inside)?;
    Some(read.seq[start_idx..end_idx].to_vec())
}
