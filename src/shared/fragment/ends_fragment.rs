use crate::{
    commands::ends::config_structs::{
        BaseQualityAggregation, BaseQualityFilter, BaseQualityFilterScope, ClipStrategy, KmerSource,
    },
    shared::{
        fragment::cigar_counts::inspect_cigar_edges,
        fragment::minimal_fragment::{
            PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
        },
        gc_tag::{GcTagValue, combine_gc_tag_values, read_gc_tag_from_record},
        indel_mode::IndelMotifFilterPolicy,
        interval::Interval,
    },
};
use anyhow::{Result, bail};
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
    /// Number of leading inside bases that still have concrete reference positions next to
    /// `boundary_pos`.
    ///
    /// This matters for `RawAlignedBoundary`: the clipped-only prefix/suffix is kept in
    /// `inside_bases`, but those bases lie outside the aligned reference span and therefore cannot
    /// be blacklist-validated against genomic coordinates.
    pub inside_reference_validation_bp: usize,
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
    pub qualities: Option<Vec<u8>>,
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
        load_base_qualities: bool,
    ) -> Result<Self> {
        let edge_info = inspect_cigar_edges(r);
        let gc_tag_value = gc_tag
            .map(|tag| read_gc_tag_from_record(r, tag))
            .unwrap_or_default();
        let qualities = if load_base_qualities {
            Some(load_record_qualities(r)?)
        } else {
            None
        };

        Ok(EndReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            left_soft_clip_bp: edge_info.left_soft_clip_bp,
            right_soft_clip_bp: edge_info.right_soft_clip_bp,
            left_motif_has_indels: motif_has_indels(
                r,
                FragmentEndSide::Left,
                clip_strategy,
                k_inside,
                edge_info.left_soft_clip_bp,
            ),
            right_motif_has_indels: motif_has_indels(
                r,
                FragmentEndSide::Right,
                clip_strategy,
                k_inside,
                edge_info.right_soft_clip_bp,
            ),
            has_hard_clip: edge_info.has_hard_clip,
            seq: r.seq().as_bytes(),
            qualities,
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
    SkipFragment,
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
    max_soft_clips: u32,
    bq_filters: &[BaseQualityFilter],
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
        ResolvedEndOutcome::SkipFragment => return None,
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
        ResolvedEndOutcome::SkipFragment => return None,
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

    let (left_end, right_end) = apply_base_quality_filters(
        forward,
        reverse,
        left_end,
        right_end,
        clip_strategy,
        k_inside,
        bq_filters,
    )?;

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
    max_soft_clips: u32,
    bq_filters: &[BaseQualityFilter],
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
        ResolvedEndOutcome::SkipFragment => return None,
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
        ResolvedEndOutcome::SkipFragment => return None,
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

    let (left_end, right_end) = apply_base_quality_filters(
        read,
        read,
        left_end,
        right_end,
        clip_strategy,
        k_inside,
        bq_filters,
    )?;

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

/// Check whether the inside-fragment motif uses aligned bases across an indel.
fn motif_has_indels(
    record: &Record,
    end_side: FragmentEndSide,
    clip_strategy: ClipStrategy,
    k_inside: usize,
    soft_clip_bp: u32,
) -> bool {
    let aligned_bases_in_motif = match clip_strategy {
        ClipStrategy::Aligned | ClipStrategy::Skip => k_inside,
        ClipStrategy::RawAlignedBoundary | ClipStrategy::RawShiftedBoundary => {
            k_inside.saturating_sub(soft_clip_bp as usize)
        }
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
/// - `ResolvedEndOutcome::SkipFragment` when the whole fragment must be discarded, e.g.
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
    max_soft_clips: u32,
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

    if soft_clip_bp > max_soft_clips {
        return ResolvedEndOutcome::SkipEndDropAssignmentBoundary;
    }

    // First decide whether indels invalidate just this end or the whole fragment.
    // Do not decide assignment-boundary behavior yet, because that depends on the clip strategy.
    if motif_has_indels && matches!(indel_filter, IndelMotifFilterPolicy::SkipAffectedFragment) {
        return ResolvedEndOutcome::SkipFragment;
    }

    let skip_end_due_to_indels = motif_has_indels
        && match indel_filter {
            IndelMotifFilterPolicy::Auto => matches!(source_inside, KmerSource::Reference),
            IndelMotifFilterPolicy::SkipAffectedEnd
            | IndelMotifFilterPolicy::SkipAffectedFragment => true,
        };

    match clip_strategy {
        ClipStrategy::Aligned | ClipStrategy::Skip => {
            if matches!(clip_strategy, ClipStrategy::Skip) && soft_clip_bp > 0 {
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
        ClipStrategy::RawAlignedBoundary | ClipStrategy::RawShiftedBoundary => {
            let assignment_boundary_pos = match end_side {
                FragmentEndSide::Left if clip_strategy.uses_shifted_boundary() => {
                    aligned_boundary_pos.saturating_sub(soft_clip_bp)
                }
                FragmentEndSide::Right if clip_strategy.uses_shifted_boundary() => {
                    aligned_boundary_pos.saturating_add(soft_clip_bp)
                }
                _ => aligned_boundary_pos,
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
    let inside_reference_validation_bp =
        inside_reference_validation_bp(read, end_side, clip_strategy, k_inside);

    Some(ResolvedFragmentEnd {
        boundary_pos,
        inside_bases,
        inside_reference_validation_bp,
    })
}

/// Apply parsed `--bq-filter` expressions to the already-resolved fragment ends.
///
/// This helper only scores the `k_inside` read bases for ends that are still
/// present after clip handling and indel filtering.
///
/// Evaluation order is intentional:
///
/// - fragment-scope filters run first on the raw candidate fragment, combining
///   the inside-base qualities from both currently resolved ends
/// - end-scope filters then run independently per surviving end and may drop
///   only the failing end
/// - if both ends are removed, the whole fragment is dropped because there is
///   nothing left to count
///
/// Returns `None` when:
///
/// - a fragment-scope filter fails, including the degenerate case where no
///   quality bases are available across all surviving ends
/// - an end-scope filter removes both ends
///
/// When `bq_filters` is empty, this is a no-op and returns the input ends
/// unchanged.
fn apply_base_quality_filters(
    left_read: &EndReadInfo,
    right_read: &EndReadInfo,
    mut left_end: Option<ResolvedFragmentEnd>,
    mut right_end: Option<ResolvedFragmentEnd>,
    clip_strategy: ClipStrategy,
    k_inside: usize,
    bq_filters: &[BaseQualityFilter],
) -> Option<(Option<ResolvedFragmentEnd>, Option<ResolvedFragmentEnd>)> {
    if bq_filters.is_empty() {
        return Some((left_end, right_end));
    }

    let left_qualities = left_end.as_ref().and_then(|_| {
        extract_inside_qualities(left_read, FragmentEndSide::Left, clip_strategy, k_inside)
    });
    let right_qualities = right_end.as_ref().and_then(|_| {
        extract_inside_qualities(right_read, FragmentEndSide::Right, clip_strategy, k_inside)
    });

    for &filter in bq_filters
        .iter()
        .filter(|filter| matches!(filter.scope, BaseQualityFilterScope::Fragment))
    {
        let fragment_score = fragment_quality_score(
            left_qualities.as_deref(),
            right_qualities.as_deref(),
            filter.aggregation,
        )?;
        if !filter.passes_value(fragment_score) {
            return None;
        }
    }

    for &filter in bq_filters
        .iter()
        .filter(|filter| matches!(filter.scope, BaseQualityFilterScope::End))
    {
        if let Some(qualities) = left_qualities.as_deref() {
            let left_score = aggregate_base_qualities(qualities, filter.aggregation)?;
            if !filter.passes_value(left_score) {
                left_end = None;
            }
        }
        if let Some(qualities) = right_qualities.as_deref() {
            let right_score = aggregate_base_qualities(qualities, filter.aggregation)?;
            if !filter.passes_value(right_score) {
                right_end = None;
            }
        }
    }

    if left_end.is_none() && right_end.is_none() {
        None
    } else {
        Some((left_end, right_end))
    }
}

fn fragment_quality_score(
    left_qualities: Option<&[u8]>,
    right_qualities: Option<&[u8]>,
    aggregation: BaseQualityAggregation,
) -> Option<f32> {
    match aggregation {
        BaseQualityAggregation::Min => left_qualities
            .iter()
            .chain(right_qualities.iter())
            .flat_map(|qualities| qualities.iter().copied())
            .min()
            .map(f32::from),
        BaseQualityAggregation::Max => left_qualities
            .iter()
            .chain(right_qualities.iter())
            .flat_map(|qualities| qualities.iter().copied())
            .max()
            .map(f32::from),
        BaseQualityAggregation::Mean => {
            let (sum, count) = left_qualities
                .iter()
                .chain(right_qualities.iter())
                .flat_map(|qualities| qualities.iter().copied())
                .fold((0_u64, 0_u64), |(sum, count), value| {
                    (sum + u64::from(value), count + 1)
                });
            (count > 0).then(|| sum as f32 / count as f32)
        }
    }
}

fn aggregate_base_qualities(qualities: &[u8], aggregation: BaseQualityAggregation) -> Option<f32> {
    if qualities.is_empty() {
        return None;
    }

    match aggregation {
        BaseQualityAggregation::Min => qualities.iter().copied().min().map(f32::from),
        BaseQualityAggregation::Mean => {
            let sum: u64 = qualities.iter().map(|&value| u64::from(value)).sum();
            Some(sum as f32 / qualities.len() as f32)
        }
        BaseQualityAggregation::Max => qualities.iter().copied().max().map(f32::from),
    }
}

/// Return how many of the extracted inside bases can still be checked against reference-backed
/// genomic positions.
///
/// Most clip strategies keep the motif boundary aligned with the genomic start/end they use for
/// validation, so all `k_inside` bases remain reference-addressable. `RawAlignedBoundary` is the
/// exception: it keeps the aligned genomic boundary but prepends/appends clipped read bases to the
/// inside motif. Those clipped-only bases do not correspond to reference positions and must be
/// excluded from blacklist validation.
fn inside_reference_validation_bp(
    read: &EndReadInfo,
    end_side: FragmentEndSide,
    clip_strategy: ClipStrategy,
    k_inside: usize,
) -> usize {
    match clip_strategy {
        ClipStrategy::Aligned | ClipStrategy::Skip | ClipStrategy::RawShiftedBoundary => k_inside,
        ClipStrategy::RawAlignedBoundary => {
            let soft_clip_bp = match end_side {
                FragmentEndSide::Left => read.left_soft_clip_bp,
                FragmentEndSide::Right => read.right_soft_clip_bp,
            };
            k_inside.saturating_sub(soft_clip_bp as usize)
        }
    }
}

fn extract_left_inside_bases(
    read: &EndReadInfo,
    clip_strategy: ClipStrategy,
    k_inside: usize,
) -> Option<Vec<u8>> {
    let (start_idx, end_idx) = inside_slice_bounds(
        read.seq.len(),
        read.left_soft_clip_bp,
        read.right_soft_clip_bp,
        FragmentEndSide::Left,
        clip_strategy,
        k_inside,
    )?;
    Some(read.seq[start_idx..end_idx].to_vec())
}

fn extract_right_inside_bases(
    read: &EndReadInfo,
    clip_strategy: ClipStrategy,
    k_inside: usize,
) -> Option<Vec<u8>> {
    let (start_idx, end_idx) = inside_slice_bounds(
        read.seq.len(),
        read.left_soft_clip_bp,
        read.right_soft_clip_bp,
        FragmentEndSide::Right,
        clip_strategy,
        k_inside,
    )?;
    Some(read.seq[start_idx..end_idx].to_vec())
}

fn extract_inside_qualities(
    read: &EndReadInfo,
    end_side: FragmentEndSide,
    clip_strategy: ClipStrategy,
    k_inside: usize,
) -> Option<Vec<u8>> {
    let qualities = read.qualities.as_ref()?;
    let (start_idx, end_idx) = inside_slice_bounds(
        qualities.len(),
        read.left_soft_clip_bp,
        read.right_soft_clip_bp,
        end_side,
        clip_strategy,
        k_inside,
    )?;
    Some(qualities[start_idx..end_idx].to_vec())
}

fn load_record_qualities(record: &Record) -> Result<Vec<u8>> {
    let qualities = record.qual().to_vec();
    if qualities.iter().all(|&value| value == 255) {
        let qname = String::from_utf8_lossy(record.qname());
        bail!(
            "BAM record '{qname}' at tid {}, pos {} has missing base qualities (`*` / 255 placeholder), but `--bq-filter` requires concrete read base qualities",
            record.tid(),
            record.pos()
        );
    }
    Ok(qualities)
}

fn inside_slice_bounds(
    len: usize,
    left_soft_clip_bp: u32,
    right_soft_clip_bp: u32,
    end_side: FragmentEndSide,
    clip_strategy: ClipStrategy,
    k_inside: usize,
) -> Option<(usize, usize)> {
    match end_side {
        FragmentEndSide::Left => {
            let start_idx = match clip_strategy {
                ClipStrategy::Aligned | ClipStrategy::Skip => left_soft_clip_bp as usize,
                ClipStrategy::RawAlignedBoundary | ClipStrategy::RawShiftedBoundary => 0,
            };
            let end_idx = start_idx.checked_add(k_inside)?;
            (end_idx <= len).then_some((start_idx, end_idx))
        }
        FragmentEndSide::Right => {
            let end_idx = match clip_strategy {
                ClipStrategy::Aligned | ClipStrategy::Skip => {
                    len.checked_sub(right_soft_clip_bp as usize)?
                }
                ClipStrategy::RawAlignedBoundary | ClipStrategy::RawShiftedBoundary => len,
            };
            let start_idx = end_idx.checked_sub(k_inside)?;
            Some((start_idx, end_idx))
        }
    }
}

#[cfg(test)]
mod tests {
    include!("ends_fragment_tests.rs");
}
