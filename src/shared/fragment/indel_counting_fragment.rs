use fxhash::FxHashMap;
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::{Cigar, Record};

use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::interval::{Interval, TouchingMergePolicy, merge_sorted_intervals};

/// Insertion anchored at one reference position with a positive inserted length.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct InsertionAnchor {
    pub reference_position: u32,
    pub inserted_length: u32,
}

/// Compact per-read info with extracted indel events.
#[derive(Debug, Clone)]
pub struct IndelReadInfo {
    pub tid: i32,
    pub interval: Interval<u32>, // Aligned reference span [start: pos(), end: reference_end())
    pub is_reverse: bool,
    /// Deletions (and ref-skips if present) as reference intervals [start, end)
    pub deletions: Vec<Interval<u32>>,
    /// Insertions anchored at one reference position with their inserted length.
    pub insertions: Vec<InsertionAnchor>,
}

impl TryFrom<&Record> for IndelReadInfo {
    type Error = crate::Error;

    #[inline]
    fn try_from(r: &Record) -> Result<Self> {
        let mut deletions: Vec<Interval<u32>> = Vec::new();
        let mut insertions: Vec<InsertionAnchor> = Vec::new();

        // Walk the CIGAR in reference coordinates
        let mut ref_pos: u32 = r.pos() as u32;

        for op in r.cigar().iter() {
            match *op {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                    ref_pos = ref_pos.saturating_add(l);
                }
                Cigar::Del(l) => {
                    let deletion_start = ref_pos;
                    let deletion_end = ref_pos.saturating_add(l);
                    if let Ok(deletion) = Interval::new(deletion_start, deletion_end) {
                        deletions.push(deletion);
                    }
                    ref_pos = deletion_end; // D consumes reference
                }
                Cigar::RefSkip(l) => {
                    // Rare in cfDNA; treat as a deletion on the reference
                    let skipped_start = ref_pos;
                    let skipped_end = ref_pos.saturating_add(l);
                    if let Ok(skipped_interval) = Interval::new(skipped_start, skipped_end) {
                        deletions.push(skipped_interval);
                    }
                    ref_pos = skipped_end; // N consumes reference
                }
                Cigar::Ins(l) => {
                    // Insertion anchored at current ref_pos
                    if l > 0 {
                        insertions.push(InsertionAnchor {
                            reference_position: ref_pos,
                            inserted_length: l,
                        });
                    }
                    // I does not consume reference
                }
                Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {
                    // Ignore: do not consume reference, and clips are not molecule
                }
            }
        }

        // Merge adjacent/overlapping deletion intervals to normalize
        if deletions.len() > 1 {
            deletions.sort_unstable_by_key(|deletion| deletion.start());
            deletions = merge_sorted_intervals(deletions, TouchingMergePolicy::MergeTouching);
        }

        Ok(IndelReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            deletions,
            insertions,
        })
    }
}

impl IndelReadInfo {
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

    /// Return the read's aligned reference span `[start, end)`.
    #[inline]
    pub fn aligned_interval(&self) -> Interval<u32> {
        self.interval
    }
}

impl PairOrientable for IndelReadInfo {
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

/// Fragment plus molecule-leaning indel accounting components.
///
/// - In non-overlap: count all per-read D/N (as deletions) and I (as insertions).
/// - In overlap: count only **supported-by-both** events:
///     * Deletions: add the intersection of deletion spans from the two reads
///     * Insertions: add at reference positions present in **both** mates (min length if disagree)
#[derive(Debug, Clone)]
pub struct FragmentWithIndelCounts {
    pub tid: i32,
    pub interval: Interval<u32>, // forward.pos .. reverse.end

    // Totals accumulated under the "pair-supported in overlap" policy:
    pub deletions_nonoverlap: u32,
    pub insertions_nonoverlap: u32,
    pub deletions_overlap_supported: u32,
    pub insertions_overlap_supported: u32,
}

impl FragmentWithIndelCounts {
    /// Inclusive fragment start on the reference.
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Exclusive fragment end on the reference.
    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Reference-span fragment length (end - start).
    #[inline]
    pub fn len_ref(&self) -> u32 {
        self.interval.len()
    }

    /// Indel-aware length: len_ref + inserts_total - deletions_total (saturating at 0).
    #[inline]
    pub fn len_indel_adjusted(&self) -> u32 {
        let ins = (self.insertions_nonoverlap as u64) + (self.insertions_overlap_supported as u64);
        let del = (self.deletions_nonoverlap as u64) + (self.deletions_overlap_supported as u64);
        let base = self.len_ref() as u64;
        base.saturating_add(ins).saturating_sub(del) as u32
    }
}

/// Partition one deletion interval into fragment-supported non-overlap bases and
/// the clipped piece that falls inside the aligned mate overlap.
fn partition_deletion_by_aligned_overlap(
    deletion_interval: Interval<u32>,
    fragment_interval: Interval<u32>,
    aligned_overlap_interval: Option<Interval<u32>>,
    nonoverlap_bases_bp: &mut u32,
    overlap_deletion_intervals: &mut Vec<Interval<u32>>,
) {
    if let Some(deletion_interval) = deletion_interval.clip_to(fragment_interval) {
        if let Some(aligned_overlap_interval) = aligned_overlap_interval {
            if let Some(left_nonoverlap_interval) =
                deletion_interval.clip_upper(aligned_overlap_interval.start())
            {
                *nonoverlap_bases_bp =
                    nonoverlap_bases_bp.saturating_add(left_nonoverlap_interval.len());
            }

            if let Some(overlap_deletion_interval) =
                deletion_interval.clip_to(aligned_overlap_interval)
            {
                overlap_deletion_intervals.push(overlap_deletion_interval);
            }

            if let Some(right_nonoverlap_interval) =
                deletion_interval.clip_lower(aligned_overlap_interval.end())
            {
                *nonoverlap_bases_bp =
                    nonoverlap_bases_bp.saturating_add(right_nonoverlap_interval.len());
            }
        } else {
            // No mate overlap at all: whole deletion is non-overlap.
            *nonoverlap_bases_bp = nonoverlap_bases_bp.saturating_add(deletion_interval.len());
        }
    }
}

/// Partition one insertion anchor into fragment-supported non-overlap bases or
/// overlap anchors keyed by reference position.
fn partition_insertion_by_aligned_overlap(
    insertion_anchor: InsertionAnchor,
    fragment_interval: Interval<u32>,
    aligned_overlap_interval: Option<Interval<u32>>,
    nonoverlap_bases_bp: &mut u32,
    overlap_insertions_by_anchor: &mut FxHashMap<u32, u32>,
) {
    let insertion_anchor_bp = insertion_anchor.reference_position;
    let inserted_length_bp = insertion_anchor.inserted_length;
    // Ignore insertions whose reference anchor lies outside the fragment span
    if !fragment_interval.contains_point(insertion_anchor_bp) {
        return;
    }
    if aligned_overlap_interval
        .is_none_or(|overlap_interval| !overlap_interval.contains_point(insertion_anchor_bp))
    {
        *nonoverlap_bases_bp = nonoverlap_bases_bp.saturating_add(inserted_length_bp);
    } else {
        // At the same ref position, keep the maximum length per read
        overlap_insertions_by_anchor
            .entry(insertion_anchor_bp)
            // Safeguards against weird cigar strings
            .and_modify(|length_bp| *length_bp = (*length_bp).max(inserted_length_bp))
            .or_insert(inserted_length_bp);
    }
}

/// Build a `FragmentWithIndelCounts` from a single read.
///
/// The fragment span is the aligned reference span of the read `[pos, reference_end)`.
/// Indel handling follows `collect_fragment_with_indel_counts` semantics:
/// - Skip when `skip_indels` is true and the read has any insertions or deletions.
/// - When `count_indels` is true, insertions increase the length and deletions decrease it.
pub fn collect_fragment_with_indel_counts_from_single_read(
    read: &IndelReadInfo,
    skip_indels: bool,
    count_indels: bool,
) -> Option<FragmentWithIndelCounts> {
    let fragment_has_indels = !read.deletions.is_empty() || !read.insertions.is_empty();
    if skip_indels && fragment_has_indels {
        return None;
    }

    let fragment_interval = read.aligned_interval();

    if !fragment_has_indels || !count_indels {
        return Some(FragmentWithIndelCounts {
            tid: read.tid,
            interval: fragment_interval,
            deletions_nonoverlap: 0,
            insertions_nonoverlap: 0,
            deletions_overlap_supported: 0,
            insertions_overlap_supported: 0,
        });
    }

    let deletions_bp: u32 = read.deletions.iter().map(|deletion| deletion.len()).sum();
    let insertions_bp: u32 = read.insertions.iter().map(|ins| ins.inserted_length).sum();

    Some(FragmentWithIndelCounts {
        tid: read.tid,
        interval: fragment_interval,
        deletions_nonoverlap: deletions_bp,
        insertions_nonoverlap: insertions_bp,
        deletions_overlap_supported: 0,
        insertions_overlap_supported: 0,
    })
}

/// Build a `FragmentWithIndelCounts` from two per-read summaries, using a
/// molecule-leaning, mate-supported policy for indel adjustments.
///
/// Concept
/// -------
/// 1) Require same contig, opposite strands, and **inward** geometry
///    (`forward.pos <= reverse.pos`). The fragment span is
///    `[forward.pos, reverse.end)` (end-exclusive).
/// 2) Split each read's indels into:
///    - **Non-overlap** (bases covered by only one mate): count fully per read
///      * Deletions/RefSkips (D/N) add to `deletions_nonoverlap`.
///      * Insertions (I)           add to `insertions_nonoverlap`.
///    - **Overlap** (bases covered by both mates): count **only if supported by both**
///      * Deletions: accumulate the **intersection** of deletion intervals
///        across the two mates -> `deletions_overlap_supported`.
///      * Insertions: count only positions where **both** mates have an insertion;
///        add `min(len_a, len_b)` at each shared reference position
///        -> `insertions_overlap_supported`.
///
/// The **indel-adjusted** length can then be derived as:
/// `len_indel_adjusted = (end - start) + insertions_total - deletions_total`,
/// where totals are the sums of non-overlap and supported-overlap components.
///
/// Parameters
/// ----------
/// - `a`, `b`: Per-read summaries with aligned reference bounds and extracted
///   indels (`IndelReadInfo`).
/// - `skip_indels`: Return `None` if a fragment has any insertions or deletions.
/// - `count_indels`: Whether to count the indels or set them to 0.
///
/// Returns
/// -------
/// - `Some(FragmentWithIndelCounts)` if the pair is inward on the same contig;
///   otherwise `None`.
pub fn collect_fragment_with_indel_counts(
    a: &IndelReadInfo,
    b: &IndelReadInfo,
    skip_indels: bool,
    count_indels: bool,
) -> Option<FragmentWithIndelCounts> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }

    // Check if fragment has any indels
    let fragment_has_indels = !forward.deletions.is_empty()
        || !forward.insertions.is_empty()
        || !reverse.deletions.is_empty()
        || !reverse.insertions.is_empty();

    if skip_indels && fragment_has_indels {
        return None;
    }

    let fragment_interval = Interval::new(forward.start(), reverse.end()).ok()?;

    // Fast path: if neither mate has any indels or we don't want to count indels,
    // return plain fragment with zero adjustments
    if !fragment_has_indels || !count_indels {
        return Some(FragmentWithIndelCounts {
            tid: forward.tid,
            interval: fragment_interval,
            deletions_nonoverlap: 0,
            insertions_nonoverlap: 0,
            deletions_overlap_supported: 0,
            insertions_overlap_supported: 0,
        });
    }

    // Reference overlap of the **aligned segments** (not the template):
    // overlap = [max(forward.start(), reverse.start()), min(forward.end(), reverse.end()))
    let aligned_overlap_start_bp = forward.start().max(reverse.start());
    let aligned_overlap_end_bp = forward.end().min(reverse.end());
    let aligned_overlap_interval =
        Interval::new(aligned_overlap_start_bp, aligned_overlap_end_bp).ok();

    // Deletions (and ref-skips)
    // Split each deletion interval into non-overlap and (possible) overlap part
    let mut deletions_nonoverlap_bp: u32 = 0;
    let mut deletion_intervals_in_overlap_forward: Vec<Interval<u32>> = Vec::new();
    let mut deletion_intervals_in_overlap_reverse: Vec<Interval<u32>> = Vec::new();

    // Extract deletions for forward read
    for &del_iv in &forward.deletions {
        partition_deletion_by_aligned_overlap(
            del_iv,
            fragment_interval,
            aligned_overlap_interval,
            &mut deletions_nonoverlap_bp,
            &mut deletion_intervals_in_overlap_forward,
        );
    }

    // Extract deletions for reverse read
    for &del_iv in &reverse.deletions {
        partition_deletion_by_aligned_overlap(
            del_iv,
            fragment_interval,
            aligned_overlap_interval,
            &mut deletions_nonoverlap_bp,
            &mut deletion_intervals_in_overlap_reverse,
        );
    }

    // Supported overlap deletions = sum of pairwise intersections
    let deletions_overlap_supported_bp = calculate_deletion_in_overlap(
        deletion_intervals_in_overlap_forward,
        deletion_intervals_in_overlap_reverse,
    );

    // Insertions
    // Non-overlap: count fully. Overlap: only if both mates have an insertion at the same ref pos.
    let mut insertions_nonoverlap_bp: u32 = 0;
    let mut insertions_in_overlap_forward: FxHashMap<u32, u32> = FxHashMap::default();
    let mut insertions_in_overlap_reverse: FxHashMap<u32, u32> = FxHashMap::default();

    // Extract insertions for forward read
    for &ins in &forward.insertions {
        partition_insertion_by_aligned_overlap(
            ins,
            fragment_interval,
            aligned_overlap_interval,
            &mut insertions_nonoverlap_bp,
            &mut insertions_in_overlap_forward,
        );
    }

    // Extract insertions for reverse read
    for &ins in &reverse.insertions {
        partition_insertion_by_aligned_overlap(
            ins,
            fragment_interval,
            aligned_overlap_interval,
            &mut insertions_nonoverlap_bp,
            &mut insertions_in_overlap_reverse,
        );
    }

    // Calculate overlap insertions
    // Both reads must agree on the position in the overlap (min insertion size of the two)
    let mut insertions_overlap_supported_bp: u32 = 0;
    if !insertions_in_overlap_forward.is_empty() && !insertions_in_overlap_reverse.is_empty() {
        for (ref_pos, len_forward) in insertions_in_overlap_forward {
            if let Some(&len_reverse) = insertions_in_overlap_reverse.get(&ref_pos) {
                insertions_overlap_supported_bp =
                    insertions_overlap_supported_bp.saturating_add(len_forward.min(len_reverse));
            }
        }
    }

    Some(FragmentWithIndelCounts {
        tid: forward.tid,
        interval: fragment_interval,
        deletions_nonoverlap: deletions_nonoverlap_bp,
        insertions_nonoverlap: insertions_nonoverlap_bp,
        deletions_overlap_supported: deletions_overlap_supported_bp,
        insertions_overlap_supported: insertions_overlap_supported_bp,
    })
}

/// Assumes intervals are start-sorted.
fn calculate_deletion_in_overlap(
    deletion_intervals_in_overlap_forward: Vec<Interval<u32>>,
    deletion_intervals_in_overlap_reverse: Vec<Interval<u32>>,
) -> u32 {
    // Supported overlap deletions.
    // Fast path for tiny lists, otherwise linear two-pointer sweep over already-sorted lists.
    let mut supported_overlap_deletions_bp: u32 = 0;

    let a = &deletion_intervals_in_overlap_forward;
    let b = &deletion_intervals_in_overlap_reverse;

    if !a.is_empty() && !b.is_empty() {
        // Tiny lists: nested loop is cheapest
        if a.len() <= 2 && b.len() <= 2 {
            for forward_deletion in a {
                for reverse_deletion in b {
                    if let Some(shared_deletion_interval) =
                        forward_deletion.intersection(*reverse_deletion)
                    {
                        supported_overlap_deletions_bp = supported_overlap_deletions_bp
                            .saturating_add(shared_deletion_interval.len());
                    }
                }
            }
        } else {
            // Larger lists: linear sweep assuming both are already start-sorted
            let (mut i, mut j) = (0usize, 0usize);
            while i < a.len() && j < b.len() {
                let forward_deletion = a[i];
                let reverse_deletion = b[j];

                if let Some(shared_deletion_interval) =
                    forward_deletion.intersection(reverse_deletion)
                {
                    supported_overlap_deletions_bp = supported_overlap_deletions_bp
                        .saturating_add(shared_deletion_interval.len());
                }

                // Advance the interval that ends first
                if forward_deletion.end() <= reverse_deletion.end() {
                    i += 1;
                } else {
                    j += 1;
                }
            }
        }
    }

    supported_overlap_deletions_bp
}

#[cfg(test)]
mod tests {
    include!("indel_counting_fragment_tests.rs");
}
