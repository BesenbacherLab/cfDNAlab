use fxhash::FxHashMap;
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::{Cigar, Record};

use crate::utils::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};

/// Compact per-read info with extracted indel events.
#[derive(Debug, Clone)]
pub struct IndelReadInfo {
    pub tid: i32,
    pub pos: u32, // Leftmost aligned reference pos
    pub end: u32, // Exclusive rightmost aligned reference end
    pub is_reverse: bool,
    /// Deletions (and ref-skips if present) as reference intervals [start, end)
    pub deletions: Vec<(u32, u32)>,
    /// Insertions as (reference position, inserted length)
    pub insertions: Vec<(u32, u32)>,
}

impl From<&Record> for IndelReadInfo {
    #[inline]
    fn from(r: &Record) -> Self {
        let mut deletions: Vec<(u32, u32)> = Vec::new();
        let mut insertions: Vec<(u32, u32)> = Vec::new();

        // Walk the CIGAR in reference coordinates
        let mut ref_pos: u32 = r.pos() as u32;

        for op in r.cigar().iter() {
            match *op {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                    ref_pos = ref_pos.saturating_add(l as u32);
                }
                Cigar::Del(l) => {
                    let s = ref_pos;
                    let e = ref_pos.saturating_add(l as u32);
                    if e > s {
                        deletions.push((s, e));
                    }
                    ref_pos = e; // D consumes reference
                }
                Cigar::RefSkip(l) => {
                    // Rare in cfDNA; treat as a deletion on the reference
                    let s = ref_pos;
                    let e = ref_pos.saturating_add(l as u32);
                    if e > s {
                        deletions.push((s, e));
                    }
                    ref_pos = e; // N consumes reference
                }
                Cigar::Ins(l) => {
                    // Insertion anchored at current ref_pos
                    if l > 0 {
                        insertions.push((ref_pos, l as u32));
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
            deletions.sort_unstable_by_key(|&(s, _)| s);
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(deletions.len());
            for (s, e) in deletions.drain(..) {
                if let Some(last) = merged.last_mut() {
                    if s <= last.1 {
                        last.1 = last.1.max(e);
                        continue;
                    }
                }
                merged.push((s, e));
            }
            deletions = merged;
        }

        IndelReadInfo {
            tid: r.tid(),
            pos: r.pos() as u32,
            end: r.reference_end() as u32,
            is_reverse: r.is_reverse(),
            deletions,
            insertions,
        }
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
        self.pos
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
    pub start: u32, // forward.pos
    pub end: u32,   // reverse.end (end-exclusive)

    // Totals accumulated under the "pair-supported in overlap" policy:
    pub deletions_nonoverlap: u32,
    pub insertions_nonoverlap: u32,
    pub deletions_overlap_supported: u32,
    pub insertions_overlap_supported: u32,
}

impl FragmentWithIndelCounts {
    /// Reference-span fragment length (end - start).
    #[inline]
    pub fn len_ref(&self) -> u32 {
        self.end - self.start
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

/// Build a `FragmentWithIndelCounts` from two `Record`s.
#[inline]
pub fn collect_fragment_with_indel_counts_from_records(
    a: &Record,
    b: &Record,
    skip_indels: bool,
    count_indels: bool,
) -> Option<FragmentWithIndelCounts> {
    let ai = IndelReadInfo::from(a);
    let bi = IndelReadInfo::from(b);
    collect_fragment_with_indel_counts(&ai, &bi, skip_indels, count_indels)
}

/// Build a `FragmentWithIndelCounts` from two per-read summaries, using a
/// molecule-leaning, mate-supported policy for indel adjustments.
///
/// Concept
/// -------
/// 1) Require same contig, opposite strands, and **inward** geometry
///    (`forward.pos <= reverse.pos`). The fragment span is
///    `[forward.pos, reverse.end)` (end-exclusive).
/// 2) Split each read’s indels into:
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

    let fragment_start_bp = forward.pos;
    let fragment_end_bp = reverse.end;

    // Fast path: if neither mate has any indels or we don't want to count indels,
    // return plain fragment with zero adjustments
    if !fragment_has_indels || !count_indels {
        return Some(FragmentWithIndelCounts {
            tid: forward.tid,
            start: fragment_start_bp,
            end: fragment_end_bp,
            deletions_nonoverlap: 0,
            insertions_nonoverlap: 0,
            deletions_overlap_supported: 0,
            insertions_overlap_supported: 0,
        });
    }

    // Reference overlap of the **aligned segments** (not the template):
    // overlap = [max(forward.pos, reverse.pos), min(forward.end, reverse.end))
    let aligned_overlap_start_bp = forward.pos.max(reverse.pos);
    let aligned_overlap_end_bp = forward.end.min(reverse.end);
    let has_aligned_overlap = aligned_overlap_end_bp > aligned_overlap_start_bp;

    // Deletions (and ref-skips)
    // Split each deletion interval into non-overlap and (possible) overlap part
    let mut deletions_nonoverlap_bp: u32 = 0;
    let mut deletion_intervals_in_overlap_forward: Vec<(u32, u32)> = Vec::new();
    let mut deletion_intervals_in_overlap_reverse: Vec<(u32, u32)> = Vec::new();

    let split_deletion_interval =
        |del_iv: (u32, u32), nonov_acc: &mut u32, ov_sink: &mut Vec<(u32, u32)>| {
            if let Some((del_start, del_end)) =
                clip_to_fragment(del_iv.0, del_iv.1, fragment_start_bp, fragment_end_bp)
            {
                if !has_aligned_overlap {
                    // Entire deletion is non-overlap if mates don't overlap at all
                    *nonov_acc = nonov_acc.saturating_add(del_end.saturating_sub(del_start));
                    return;
                }

                // Left non-overlap segment
                if del_start < aligned_overlap_start_bp {
                    *nonov_acc = nonov_acc
                        .saturating_add(aligned_overlap_start_bp.saturating_sub(del_start));
                }

                // Overlap segment
                let ov_s = del_start.max(aligned_overlap_start_bp);
                let ov_e = del_end.min(aligned_overlap_end_bp);
                if ov_e > ov_s {
                    ov_sink.push((ov_s, ov_e));
                }

                // Right non-overlap segment
                if del_end > aligned_overlap_end_bp {
                    *nonov_acc =
                        nonov_acc.saturating_add(del_end.saturating_sub(aligned_overlap_end_bp));
                }
            }
        };

    // Extract deletions for forward read
    for &del_iv in &forward.deletions {
        split_deletion_interval(
            del_iv,
            &mut deletions_nonoverlap_bp,
            &mut deletion_intervals_in_overlap_forward,
        );
    }

    // Extract deletions for reverse read
    for &del_iv in &reverse.deletions {
        split_deletion_interval(
            del_iv,
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

    let split_insertion = |(ins_ref_pos, ins_len): (u32, u32),
                           nonov_acc: &mut u32,
                           ov_map: &mut FxHashMap<u32, u32>| {
        if !has_aligned_overlap
            || ins_ref_pos < aligned_overlap_start_bp
            || ins_ref_pos >= aligned_overlap_end_bp
        {
            *nonov_acc = nonov_acc.saturating_add(ins_len);
        } else {
            // At the same ref position, keep the maximum length per read
            ov_map
                .entry(ins_ref_pos)
                // Safeguards against weird cigar strings
                .and_modify(|x| *x = (*x).max(ins_len))
                .or_insert(ins_len);
        }
    };

    // Extract insertions for forward read
    for &ins in &forward.insertions {
        split_insertion(
            ins,
            &mut insertions_nonoverlap_bp,
            &mut insertions_in_overlap_forward,
        );
    }

    // Extract insertions for reverse read
    for &ins in &reverse.insertions {
        split_insertion(
            ins,
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
        start: fragment_start_bp,
        end: fragment_end_bp,
        deletions_nonoverlap: deletions_nonoverlap_bp,
        insertions_nonoverlap: insertions_nonoverlap_bp,
        deletions_overlap_supported: deletions_overlap_supported_bp,
        insertions_overlap_supported: insertions_overlap_supported_bp,
    })
}

/// Assumes intervals are start-sorted.
fn calculate_deletion_in_overlap(
    deletion_intervals_in_overlap_forward: Vec<(u32, u32)>,
    deletion_intervals_in_overlap_reverse: Vec<(u32, u32)>,
) -> u32 {
    // Supported overlap deletions.
    // Fast path for tiny lists, otherwise linear two-pointer sweep over already-sorted lists.
    let mut supported_overlap_deletions_bp: u32 = 0;

    let a = &deletion_intervals_in_overlap_forward;
    let b = &deletion_intervals_in_overlap_reverse;

    if !a.is_empty() && !b.is_empty() {
        // Tiny lists: nested loop is cheapest
        if a.len() <= 2 && b.len() <= 2 {
            for &(f_start, f_end) in a {
                for &(r_start, r_end) in b {
                    let start = f_start.max(r_start);
                    let end = f_end.min(r_end);
                    if end > start {
                        supported_overlap_deletions_bp =
                            supported_overlap_deletions_bp.saturating_add(end - start);
                    }
                }
            }
        } else {
            // Larger lists: linear sweep assuming both are already start-sorted
            let (mut i, mut j) = (0usize, 0usize);
            while i < a.len() && j < b.len() {
                let (f_start, f_end) = a[i];
                let (r_start, r_end) = b[j];

                let start = f_start.max(r_start);
                let end = f_end.min(r_end);
                if end > start {
                    supported_overlap_deletions_bp =
                        supported_overlap_deletions_bp.saturating_add(end - start);
                }

                // Advance the interval that ends first
                if f_end <= r_end {
                    i += 1;
                } else {
                    j += 1;
                }
            }
        }
    }

    supported_overlap_deletions_bp
}

// Clip [s,e) to the fragment span; return None if outside
#[inline]
fn clip_to_fragment(s: u32, e: u32, frag_s: u32, frag_e: u32) -> Option<(u32, u32)> {
    let cs = s.max(frag_s);
    let ce = e.min(frag_e);
    (ce > cs).then_some((cs, ce))
}
