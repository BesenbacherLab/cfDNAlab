use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::interval::Interval;
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::Record;

/// Compact per-read info with extracted indel events.
#[derive(Debug, Clone)]
pub(crate) struct FragReadInfo {
    pub(crate) tid: i32,
    pub(crate) interval: Interval<u32>, // Aligned reference span [start: pos(), end: reference_end())
    pub(crate) is_reverse: bool,
    pub(crate) mapq: u8,
    /// Aligned strand
    pub(crate) strand: char,
    /// Whether this read is read1 (or `false` for read2).
    pub(crate) is_read_1: bool,
}

impl TryFrom<&Record> for FragReadInfo {
    type Error = crate::Error;

    #[inline]
    fn try_from(r: &Record) -> Result<Self> {
        let strand: char = r.strand().strand_symbol().chars().collect::<Vec<char>>()[0];
        Ok(FragReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            mapq: r.mapq(),
            strand,
            is_read_1: r.is_first_in_template(),
        })
    }
}

impl FragReadInfo {
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub(crate) fn end(&self) -> u32 {
        self.interval.end()
    }
}

impl PairOrientable for FragReadInfo {
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

/// Fragment with mapq and read1 strand.
#[derive(Debug, Clone)]
pub(crate) struct FragFileFragment {
    pub(crate) interval: Interval<u32>, // fragment span [forward.pos, reverse.reference_end)
    pub(crate) min_mapq: u8,
    pub(crate) read1_strand: char,
}

impl FragFileFragment {
    /// Inclusive fragment start on the reference.
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Exclusive fragment end on the reference.
    #[inline]
    pub(crate) fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Reference-span fragment length (end - start).
    #[inline]
    pub(crate) fn len(&self) -> u32 {
        self.interval.len()
    }
}

/// Build a `FragFileFragment` from two reads.
pub(crate) fn collect_fragment_with_frag_file_info(
    a: &FragReadInfo,
    b: &FragReadInfo,
) -> Option<FragFileFragment> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }

    // If both reads have the same "is_read_1" status, we skip the fragment
    if !forward.is_read_1 && !reverse.is_read_1 {
        return None;
    }
    if forward.is_read_1 && reverse.is_read_1 {
        return None;
    }

    Some(FragFileFragment {
        interval: Interval::new(forward.start(), reverse.end()).ok()?,
        min_mapq: forward.mapq.min(reverse.mapq),
        read1_strand: if forward.is_read_1 {
            forward.strand
        } else {
            reverse.strand
        },
    })
}

/// Build a `FragFileFragment` from a single read (unpaired input).
pub(crate) fn collect_fragment_with_frag_file_info_from_single_read(
    read: &FragReadInfo,
) -> Option<FragFileFragment> {
    Some(FragFileFragment {
        interval: read.interval,
        min_mapq: read.mapq,
        // No read 1/2 so just the strand of the single read
        read1_strand: read.strand,
    })
}
