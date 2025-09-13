use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::Record;

/// Basic fragment on the reference (0-based, end-exclusive).
#[derive(Debug, Clone, Copy)]
pub struct Fragment {
    /// 0-based tid/contig id
    pub tid: i32,
    /// 0-based inclusive start (left boundary of the forward read)
    pub start: u32,
    /// 0-based exclusive end (right boundary of the reverse read)
    pub end: u32,
}

impl Fragment {
    /// Length of the fragment (end - start).
    pub fn len(&self) -> u32 {
        (self.end - self.start) as u32
    }
}

/// Minimal per-read info needed to build a Fragment without stashing full Records.
#[derive(Debug, Clone, Copy)]
pub struct MinimalReadInfo {
    pub tid: i32, // Contig id
    pub pos: u32, // 0-based leftmost aligned ref pos
    pub end: u32, // 0-based exclusive rightmost aligned ref pos (reference_end)
    pub is_reverse: bool,
}

impl From<&Record> for MinimalReadInfo {
    #[inline]
    fn from(r: &Record) -> Self {
        MinimalReadInfo {
            tid: r.tid(),
            pos: r.pos() as u32,
            end: r.reference_end() as u32,
            is_reverse: r.is_reverse(),
        }
    }
}

impl PairOrientable for MinimalReadInfo {
    #[inline]
    fn tid(&self) -> i32 {
        self.tid
    }
    #[inline]
    fn is_reverse(&self) -> bool {
        self.is_reverse
    }
}

/// Compute the cfDNA fragment coordinates (forward.left -> reverse.right).
///
/// Parameters
/// ----------
/// a: &Record
///     One read of the pair (mapped).
/// b: &Record
///     The mate read (mapped).
///
/// Returns
/// -------
/// frag: Option<Fragment>
///     The fragment if both reads are mapped to the same contig, on opposite strands,
///     and inward-facing; otherwise `None`.
pub fn collect_fragment_from_records(a: &Record, b: &Record) -> Option<Fragment> {
    collect_fragment(&MinimalReadInfo::from(a), &MinimalReadInfo::from(b))
}

/// Build a Fragment from two `MinimalReadInfo`s (no full BAM records needed).
pub fn collect_fragment(a: &MinimalReadInfo, b: &MinimalReadInfo) -> Option<Fragment> {
    let (fwd, rev) = oriented_pair_from_read_info(a, b)?;
    if rev.end <= fwd.pos {
        return None;
    }
    Some(Fragment {
        tid: fwd.tid,
        start: fwd.pos,
        end: rev.end,
    })
}

/* --- Helpers --- */

/// Pair-orientation trait so we can write a single generic function for orienting pairs
pub trait PairOrientable {
    fn tid(&self) -> i32;
    fn is_reverse(&self) -> bool;
}

// /// Identify forward/reverse reads (return (forward, reverse)) if both are inward.
// ///
// /// Parameters
// /// ----------
// /// a: &Record
// ///     One read.
// /// b: &Record
// ///     Mate read.
// ///
// /// Returns
// /// -------
// /// pair: Option<(&Record, &Record)>
// ///     `(forward, reverse)` or `None` if invalid (different contigs, same strand).
// fn oriented_pair_from_records<'a>(
//     a: &'a Record,
//     b: &'a Record,
// ) -> Option<(&'a Record, &'a Record)> {
//     if a.tid() != b.tid() {
//         return None;
//     }
//     match (a.is_reverse(), b.is_reverse()) {
//         (false, true) => Some((a, b)),
//         (true, false) => Some((b, a)),
//         _ => None,
//     }
// }

/// Identify forward/reverse reads (generic to PairOrientable)
/// (return (forward, reverse)) if both are inward.
///
/// Parameters
/// ----------
/// a: &MinimalReadInfo
///     One read.
/// b: &MinimalReadInfo
///     Mate read.
///
/// Returns
/// -------
/// pair: Option<(&MinimalReadInfo, &MinimalReadInfo)>
///     `(forward, reverse)` or `None` if invalid (different contigs, same strand).
#[inline]
pub fn oriented_pair_from_read_info<'a, T: PairOrientable>(
    a: &'a T,
    b: &'a T,
) -> Option<(&'a T, &'a T)> {
    if a.tid() != b.tid() {
        return None;
    }
    match (a.is_reverse(), b.is_reverse()) {
        (false, true) => Some((a, b)), // a forward, b reverse
        (true, false) => Some((b, a)), // b forward, a reverse
        _ => None,                     // same orientation or ambiguous
    }
}

// Other ideas but commented out for now!

// /// Reference-overlap between mates (0-based, end-exclusive).
// #[derive(Debug, Clone, Copy)]
// pub struct Overlap {
//     pub start: u32,
//     pub end: u32,
// }

// impl Overlap {
//     pub fn len(&self) -> u32 {
//         self.end - self.start
//     }
// }

/* --- Simple fragment --- */

/* --- Overlap (match-/mismatch-only ─── */

// /// Overlap data restricted to **matches/mismatches only**:
// /// stores only columns where *both* reads have an aligned base (M/=/X) at the same
// /// reference coordinate inside the overlap window. Insertions are **dropped**.
// ///
// /// Designed so you can do:
// /// - `left[i] == right[i]`              (pairwise equality in overlap)
// /// - `left[i] == ref[ref_coords[i]]`    (compare to reference genome if needed)
// #[derive(Debug, Clone)]
// pub struct FragmentOverlapMM {
//     pub frag: Fragment,
//     pub overlap: Overlap,
//     /// 0-based reference coordinates for each column (strictly increasing).
//     pub ref_coords: Vec<u32>,
//     /// Bases from the forward read at those coords (left -> right on reference).
//     pub left_bases: Vec<u8>,
//     /// Bases from the reverse read at those coords (left -> right on reference).
//     pub right_bases: Vec<u8>,
// }

// impl FragmentOverlapMM {
//     /// Build the MM-only overlap from a read pair.
//     ///
//     /// Parameters
//     /// ----------
//     /// a: &Record
//     ///     One read of the pair (mapped).
//     /// b: &Record
//     ///     Mate read (mapped).
//     ///
//     /// Returns
//     /// -------
//     /// ovl: Option<FragmentOverlapMM>
//     ///     `Some` with equal-length `ref_coords`, `left`, `right` if overlap exists; `None` otherwise.
//     pub fn from_pair(a: &Record, b: &Record) -> Option<Self> {
//         let (fwd, rev) = oriented_pair_from_records(a, b)?;
//         let frag = collect_fragment_from_records(fwd, rev)?;
//         let ov = overlap_of_mates(fwd, rev)?;
//         if ov.len() <= 0 {
//             return None;
//         }

//         // Extract only aligned (M/=/X) bases in the overlap window from each read.
//         let (l_bases, l_coords) = extract_aligned_in_range(fwd, (ov.start, ov.end));
//         let (r_bases, r_coords) = extract_aligned_in_range(rev, (ov.start, ov.end));

//         // Merge-join by reference coordinate to keep only coords present in both.
//         let (ref_coords, left_bases, right_bases) =
//             merge_by_coord_mm(&l_coords, &l_bases, &r_coords, &r_bases);

//         if ref_coords.is_empty() {
//             None
//         } else {
//             Some(Self {
//                 frag,
//                 overlap: ov,
//                 ref_coords,
//                 left_bases,
//                 right_bases,
//             })
//         }
//     }
// }

/* --- Sequences within-fragment --- */

// TODO: Check from cigar stats whether we need to check the cigar string?

// /// Minimal flags summarizing which CIGAR categories were present **within the fragment window**.
// #[derive(Debug, Clone, Copy, Default)]
// pub struct ReadSliceInfo {
//     pub has_insertion: bool, // I
//     pub has_deletion: bool,  // D
//     pub has_refskip: bool,   // N
//     pub has_softclip: bool,  // S (at edges adjoining the window)
// }

// /// Fragment plus per-read sequences **within the fragment** (left->right in reference orientation).
// /// Sequences include aligned bases and insertions that occur inside the fragment window.
// /// Soft-clipped bases are excluded.
// ///
// /// Use this when you need fragment-end motifs or per-read sequence context across the whole fragment.
// #[derive(Debug, Clone)]
// pub struct FragmentWithSequences {
//     pub frag: Fragment,
//     pub left_seq: Vec<u8>, // forward read within fragment (aligned + insertions)
//     pub right_seq: Vec<u8>, // reverse read within fragment (aligned + insertions)
//     pub left_info: ReadSliceInfo,
//     pub right_info: ReadSliceInfo,
// }

// impl FragmentWithSequences {
//     /// Build `FragmentWithSequences` from a read pair.
//     ///
//     /// Parameters
//     /// ----------
//     /// a: &Record
//     ///     One read of the pair (mapped).
//     /// b: &Record
//     ///     Mate read (mapped).
//     ///
//     /// Returns
//     /// -------
//     /// frag_seq: Option<FragmentWithSequences>
//     ///     `Some` if a valid inward-facing fragment exists; `None` otherwise.
//     pub fn from_pair(a: &Record, b: &Record) -> Option<Self> {
//         let (fwd, rev) = oriented_pair_from_records(a, b)?;
//         let frag = collect_fragment_from_records(fwd, rev)?;

//         // Trim each read to the part that lies within the fragment on the reference.
//         let f_range = (
//             fwd.pos() as u32,
//             fwd.reference_end().min(frag.end as i64) as u32,
//         );
//         let r_range = (
//             rev.pos().max(frag.start as i64) as u32,
//             rev.reference_end() as u32,
//         );

//         let (left_seq, left_info) = slice_read_to_range(fwd, f_range);
//         let (right_seq, right_info) = slice_read_to_range(rev, r_range);

//         Some(Self {
//             frag,
//             left_seq,
//             right_seq,
//             left_info,
//             right_info,
//         })
//     }
// }

// /// Intersection of the two read alignments on the reference (0-based, end-exclusive).
// fn overlap_of_mates(fwd: &Record, rev: &Record) -> Option<Overlap> {
//     let f_start = fwd.pos();
//     let f_end = fwd.reference_end();
//     let r_start = rev.pos();
//     let r_end = rev.reference_end();
//     let start = f_start.max(r_start);
//     let end = f_end.min(r_end);
//     if end > start {
//         Some(Overlap {
//             start: start as u32,
//             end: end as u32,
//         })
//     } else {
//         None
//     }
// }

// /// Extract only aligned (M/=/X) bases from `range` and their reference coordinates.
// ///
// /// Parameters
// /// ----------
// /// rec: &Record
// ///     BAM record (mapped).
// /// range: (u32, u32)
// ///     0-based inclusive start, exclusive end on the reference.
// ///
// /// Returns
// /// -------
// /// result: (Vec<u8>, Vec<u32>)
// ///     (bases, ref_coords) for aligned columns inside `range`.
// fn extract_aligned_in_range(rec: &Record, range: (u32, u32)) -> (Vec<u8>, Vec<u32>) {
//     let (range_start, range_end) = range;
//     if range_end <= range_start {
//         return (Vec::new(), Vec::new());
//     }

//     let seq = rec.seq().as_bytes();
//     let mut out_seq: Vec<u8> = Vec::new();
//     let mut out_pos: Vec<u32> = Vec::new();

//     let mut ref_pos = rec.pos() as u32;
//     let mut read_pos: usize = 0;

//     for op in rec.cigar().iter() {
//         match *op {
//             Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
//                 let op_ref_start = ref_pos;
//                 let op_ref_end = ref_pos + l;

//                 // Overlap with requested range on the reference
//                 let ov_start = op_ref_start.max(range_start);
//                 let ov_end = op_ref_end.min(range_end);
//                 if ov_end > ov_start {
//                     let offset = (ov_start - op_ref_start) as usize;
//                     let len = (ov_end - ov_start) as usize;
//                     let s = &seq[read_pos + offset..read_pos + offset + len];
//                     out_seq.extend_from_slice(s);
//                     for k in 0..len {
//                         out_pos.push(ov_start + k as u32);
//                     }
//                 }
//                 ref_pos += l;
//                 read_pos += l as usize;
//             }
//             Cigar::Ins(l) => {
//                 // Drop insertions entirely for MM-only paths.
//                 read_pos += l as usize;
//             }
//             Cigar::Del(l) | Cigar::RefSkip(l) => {
//                 ref_pos += l;
//             }
//             Cigar::SoftClip(l) => {
//                 read_pos += l as usize;
//             }
//             Cigar::HardClip(_) | Cigar::Pad(_) => { /* no-op */ }
//         }
//     }

//     (out_seq, out_pos)
// }

// /// Merge-join two (coords, bases) arrays into shared-coordinate columns.
// /// All outputs have equal length.
// ///
// /// Parameters
// /// ----------
// /// lc, lb, rc, rb:
// ///     Left/right reference coords and respective sequence bases, each sorted by coord.
// ///
// /// Returns
// /// -------
// /// result: (Vec<u32>, Vec<u8>, Vec<u8>)
// ///     (ref_coords, left_bases, right_bases) at positions covered by both reads.
// fn merge_by_coord_mm(lc: &[u32], lb: &[u8], rc: &[u32], rb: &[u8]) -> (Vec<u32>, Vec<u8>, Vec<u8>) {
//     let mut i = 0usize;
//     let mut j = 0usize;
//     let mut coords = Vec::new();
//     let mut left = Vec::new();
//     let mut right = Vec::new();

//     while i < lc.len() && j < rc.len() {
//         if lc[i] == rc[j] {
//             coords.push(lc[i]);
//             left.push(lb[i]);
//             right.push(rb[j]);
//             i += 1;
//             j += 1;
//         } else if lc[i] < rc[j] {
//             i += 1;
//         } else {
//             j += 1;
//         }
//     }

//     (coords, left, right)
// }

// /// Slice a read to a reference range, including aligned bases and insertions within the window.
// /// Soft-clipped bases are excluded. Also returns a lightweight summary of CIGAR categories
// /// encountered inside the window.
// ///
// /// Parameters
// /// ----------
// /// rec: &Record
// ///     BAM record (mapped).
// /// range: (u32, u32)
// ///     0-based inclusive start, exclusive end on the reference.
// ///
// /// Returns
// /// -------
// /// result: (Vec<u8>, ReadSliceInfo)
// ///     (sequence, summary flags) within the requested reference range.
// fn slice_read_to_range(rec: &Record, range: (u32, u32)) -> (Vec<u8>, ReadSliceInfo) {
//     let (range_start, range_end) = range;
//     if range_end <= range_start {
//         return (Vec::new(), ReadSliceInfo::default());
//     }

//     let seq = rec.seq().as_bytes();
//     let mut out: Vec<u8> = Vec::new();
//     let mut info = ReadSliceInfo::default();

//     let mut ref_pos = rec.pos() as u32;
//     let mut read_pos: usize = 0;

//     for op in rec.cigar().iter() {
//         match *op {
//             Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
//                 let op_ref_start = ref_pos;
//                 let op_ref_end = ref_pos + l;
//                 let ov_start = op_ref_start.max(range_start);
//                 let ov_end = op_ref_end.min(range_end);
//                 if ov_end > ov_start {
//                     let offset = (ov_start - op_ref_start) as usize;
//                     let len = (ov_end - ov_start) as usize;
//                     out.extend_from_slice(&seq[read_pos + offset..read_pos + offset + len]);
//                 }
//                 ref_pos += l;
//                 read_pos += l as usize;
//             }
//             Cigar::Ins(l) => {
//                 info.has_insertion = true;
//                 if ref_pos >= range_start && ref_pos < range_end {
//                     let len = l as usize;
//                     out.extend_from_slice(&seq[read_pos..read_pos + len]);
//                 }
//                 read_pos += l as usize;
//             }
//             Cigar::Del(l) => {
//                 info.has_deletion = true;
//                 ref_pos += l;
//             }
//             Cigar::RefSkip(l) => {
//                 info.has_refskip = true;
//                 ref_pos += l;
//             }
//             Cigar::SoftClip(l) => {
//                 info.has_softclip = true;
//                 read_pos += l as usize;
//             }
//             Cigar::HardClip(_) | Cigar::Pad(_) => { /* no-op */ }
//         }
//     }

//     (out, info)
// }
