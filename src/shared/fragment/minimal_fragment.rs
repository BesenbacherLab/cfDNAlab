use crate::Result;
use crate::shared::gc_tag::{GcTagValue, combine_gc_tag_values};
use crate::shared::interval::Interval;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Record;

/// Basic fragment on the reference (0-based, end-exclusive).
#[derive(Debug, Clone, Copy)]
pub struct Fragment {
    /// tid/contig id
    pub tid: i32,
    /// Checked non-empty fragment span on the reference.
    pub interval: Interval<u32>,
    /// Optional GC weight from aux tag if provided
    pub gc_tag: crate::shared::gc_tag::GcTagValue,
}

impl Fragment {
    /// Inclusive start (left boundary of the forward read).
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Exclusive end (right boundary of the reverse read).
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

/// Minimal per-read info needed to build a Fragment without stashing full Records.
#[derive(Debug, Clone, Copy)]
pub struct MinimalReadInfo {
    pub tid: i32,                // Contig id
    pub interval: Interval<u32>, // Aligned reference span [start: pos(), end: reference_end())
    pub is_reverse: bool,
    pub gc_tag: crate::shared::gc_tag::GcTagValue,
}

impl TryFrom<&Record> for MinimalReadInfo {
    type Error = crate::Error;

    #[inline]
    fn try_from(r: &Record) -> Result<Self> {
        Ok(MinimalReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            // Default tag value so tag reads are opt-in at iterator level
            gc_tag: GcTagValue::default(),
        })
    }
}

impl MinimalReadInfo {
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
        self.start()
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
    collect_fragment(
        &MinimalReadInfo::try_from(a).ok()?,
        &MinimalReadInfo::try_from(b).ok()?,
    )
}

/// Build a Fragment from two `MinimalReadInfo`s (no full BAM records needed).
pub fn collect_fragment(a: &MinimalReadInfo, b: &MinimalReadInfo) -> Option<Fragment> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }
    let gc_tag = combine_gc_tag_values(&forward.gc_tag, &reverse.gc_tag);
    let fragment_interval = Interval::new(forward.start(), reverse.end()).ok()?;
    Some(Fragment {
        tid: forward.tid,
        interval: fragment_interval,
        gc_tag,
    })
}

/// Build a Fragment from a single read (unpaired input).
pub fn collect_fragment_from_single_read(read: &MinimalReadInfo) -> Option<Fragment> {
    Some(Fragment {
        tid: read.tid,
        interval: read.aligned_interval(),
        gc_tag: read.gc_tag,
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
