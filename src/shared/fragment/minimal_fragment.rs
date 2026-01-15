use crate::shared::gc_tag::{GcTagValue, combine_gc_tag_values};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Record;

/// Basic fragment on the reference (0-based, end-exclusive).
#[derive(Debug, Clone, Copy)]
pub struct Fragment {
    /// tid/contig id
    pub tid: i32,
    /// inclusive start (left boundary of the forward read)
    pub start: u32,
    /// exclusive end (right boundary of the reverse read)
    pub end: u32,
    /// Optional GC weight from aux tag if provided
    pub gc_tag: crate::shared::gc_tag::GcTagValue,
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
    pub pos: u32, // leftmost aligned ref pos
    pub end: u32, // exclusive rightmost aligned ref pos (reference_end)
    pub is_reverse: bool,
    pub gc_tag: crate::shared::gc_tag::GcTagValue,
}

impl From<&Record> for MinimalReadInfo {
    #[inline]
    fn from(r: &Record) -> Self {
        MinimalReadInfo {
            tid: r.tid(),
            pos: r.pos() as u32,
            end: r.reference_end() as u32,
            is_reverse: r.is_reverse(),
            gc_tag: GcTagValue::default(),
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
    #[inline]
    fn pos(&self) -> u32 {
        self.pos
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
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }
    let gc_tag = combine_gc_tag_values(&forward.gc_tag, &reverse.gc_tag);
    Some(Fragment {
        tid: forward.tid,
        start: forward.pos,
        end: reverse.end,
        gc_tag,
    })
}

/* --- Helpers --- */

/// Pair-orientation trait so we can write a single generic function for orienting pairs
pub trait PairOrientable {
    fn tid(&self) -> i32;
    fn is_reverse(&self) -> bool;
    fn pos(&self) -> u32;
}

/// Identify forward/reverse reads (generic to PairOrientable)
/// (return (forward, reverse)) if both are inward.
///
/// Parameters
/// ----------
/// a:
///     One read.
/// b:
///     Mate read.
///
/// Returns
/// -------
/// pair: `(forward, reverse)` or `None` if invalid (different contigs, same strand).
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

/// Whether a fragment from two reads are inwards-oriented, meaning `forward.pos <= reverse.pos`.
///
/// Parameters
/// ----------
/// forward:
///     The forward read.
/// reverse:
///     The reverse read.
#[inline]
pub fn is_inwards_oriented<'a, T: PairOrientable>(forward: &'a T, backward: &'a T) -> bool {
    forward.pos() <= backward.pos()
}
