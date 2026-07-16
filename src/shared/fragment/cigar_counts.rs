use rust_htslib::bam::{Record, record::Cigar};

#[cfg(feature = "cmd_lengths")]
use crate::shared::interval::{Interval, TouchingMergePolicy, merge_sorted_intervals};

/// Terminal clipping summary in BAM storage orientation.
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub(crate) struct CigarEdgeInfo {
    pub(crate) left_soft_clip_bp: u32,
    pub(crate) right_soft_clip_bp: u32,
    pub(crate) has_hard_clip: bool,
}

/// Insertion anchored at one reference position with a positive inserted length.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
#[cfg(feature = "cmd_lengths")]
pub(crate) struct InsertionAnchor {
    pub(crate) reference_position: u32,
    pub(crate) inserted_length: u32,
}

/// Compact indel summary extracted from one record's CIGAR string.
#[derive(Debug, Clone, Default, Eq, PartialEq)]
#[cfg(feature = "cmd_lengths")]
pub(crate) struct CigarIndelInfo {
    /// Deletions and ref-skips as reference intervals `[start, end)`.
    pub(crate) deletions: Vec<Interval<u32>>,
    /// Insertions anchored at one reference position with their inserted length.
    pub(crate) insertions: Vec<InsertionAnchor>,
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
pub(crate) fn inspect_cigar_edges(record: &Record) -> CigarEdgeInfo {
    let mut info = CigarEdgeInfo::default();

    for op in record.cigar().iter() {
        match *op {
            Cigar::SoftClip(bp) => {
                info.left_soft_clip_bp = info.left_soft_clip_bp.saturating_add(bp);
            }
            Cigar::HardClip(_) => {
                info.has_hard_clip = true;
            }
            _ => break,
        }
    }

    for op in record.cigar().iter().rev() {
        match *op {
            Cigar::SoftClip(bp) => {
                info.right_soft_clip_bp = info.right_soft_clip_bp.saturating_add(bp);
            }
            Cigar::HardClip(_) => {
                info.has_hard_clip = true;
            }
            _ => break,
        }
    }

    info
}

/// Extract deletions/ref-skips and insertions from a record's CIGAR string.
///
/// Reference-consuming operations are tracked in aligned reference coordinates starting at
/// `record.pos()`. Adjacent or overlapping deletion-like intervals are merged after parsing so
/// downstream code can assume a normalized representation.
#[cfg(feature = "cmd_lengths")]
pub(crate) fn inspect_cigar_indels(record: &Record) -> CigarIndelInfo {
    let mut deletions: Vec<Interval<u32>> = Vec::new();
    let mut insertions: Vec<InsertionAnchor> = Vec::new();

    // Walk the CIGAR in reference coordinates
    let mut ref_pos: u32 = record.pos() as u32;

    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(length_bp) | Cigar::Equal(length_bp) | Cigar::Diff(length_bp) => {
                ref_pos = ref_pos.saturating_add(length_bp);
            }
            Cigar::Del(length_bp) => {
                let deletion_start = ref_pos;
                let deletion_end = ref_pos.saturating_add(length_bp);
                if let Ok(deletion) = Interval::new(deletion_start, deletion_end) {
                    deletions.push(deletion);
                }
                ref_pos = deletion_end; // D consumes reference
            }
            Cigar::RefSkip(length_bp) => {
                // Rare in cfDNA; treat as a deletion on the reference.
                let skipped_start = ref_pos;
                let skipped_end = ref_pos.saturating_add(length_bp);
                if let Ok(skipped_interval) = Interval::new(skipped_start, skipped_end) {
                    deletions.push(skipped_interval);
                }
                ref_pos = skipped_end; // N consumes reference
            }
            Cigar::Ins(length_bp) => {
                // Insertion anchored at current ref_pos.
                if length_bp > 0 {
                    insertions.push(InsertionAnchor {
                        reference_position: ref_pos,
                        inserted_length: length_bp,
                    });
                }
                // I does not consume reference
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {
                // Ignore: do not consume reference, and clips are not molecule
            }
        }
    }

    // Merge adjacent/overlapping deletion intervals to normalize.
    if deletions.len() > 1 {
        deletions.sort_unstable_by_key(|deletion| deletion.start());
        deletions = merge_sorted_intervals(deletions, TouchingMergePolicy::MergeTouching);
    }

    CigarIndelInfo {
        deletions,
        insertions,
    }
}

#[cfg(test)]
mod tests {
    include!("cigar_counts_tests.rs");
}
