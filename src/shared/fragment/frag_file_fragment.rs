use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::interval::Interval;
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::Record;

/// Compact per-read info with extracted indel events.
#[derive(Debug, Clone)]
pub struct FragReadInfo {
    pub tid: i32,
    pub interval: Interval<u32>, // Aligned reference span [start: pos(), end: reference_end())
    pub is_reverse: bool,
    pub mapq: u8,
    /// Aligned strand
    pub strand: char,
    /// Whether this read is read1 (or `false` for read2).
    pub is_read_1: bool,
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
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub fn end(&self) -> u32 {
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
pub struct FragFileFragment {
    pub tid: i32,
    pub interval: Interval<u32>, // fragment span [forward.pos, reverse.reference_end)
    pub min_mapq: u8,
    pub read1_strand: char,
}

impl FragFileFragment {
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
    pub fn len(&self) -> u32 {
        self.interval.len()
    }
}

/// Build a `FragFileFragment` from two reads.
pub fn collect_fragment_with_frag_file_info(
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
        tid: forward.tid,
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
pub fn collect_fragment_with_frag_file_info_from_single_read(
    read: &FragReadInfo,
) -> Option<FragFileFragment> {
    Some(FragFileFragment {
        tid: read.tid,
        interval: read.interval,
        min_mapq: read.mapq,
        // No read 1/2 so just the strand of the single read
        read1_strand: read.strand,
    })
}
