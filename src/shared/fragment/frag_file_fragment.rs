use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::Record;

/// Compact per-read info with extracted indel events.
#[derive(Debug, Clone)]
pub struct FragReadInfo {
    pub tid: i32,
    pub pos: u32, // Leftmost aligned reference pos
    pub end: u32, // Exclusive rightmost aligned reference end
    pub is_reverse: bool,
    pub mapq: u8,
    /// Aligned strand
    pub strand: char,
    /// Whether this read is read1 (or `false` for read2).
    pub is_read_1: bool,
}

impl From<&Record> for FragReadInfo {
    #[inline]
    fn from(r: &Record) -> Self {
        let strand: char = r.strand().strand_symbol().chars().collect::<Vec<char>>()[0];
        FragReadInfo {
            tid: r.tid(),
            pos: r.pos() as u32,
            end: r.reference_end() as u32,
            is_reverse: r.is_reverse(),
            mapq: r.mapq(),
            strand: strand,
            is_read_1: r.is_first_in_template(),
        }
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
        self.pos
    }
}

/// Fragment with mapq and read1 strand.
#[derive(Debug, Clone)]
pub struct FragFileFragment {
    pub tid: i32,
    pub start: u32, // forward.pos
    pub end: u32,   // reverse.end (end-exclusive)
    pub min_mapq: u8,
    pub read1_strand: char,
}

impl FragFileFragment {
    /// Reference-span fragment length (end - start).
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start
    }
}

/// Build a `FragFileFragment` from two `Record`s.
#[inline]
pub fn collect_fragment_with_frag_file_info_from_records(
    a: &Record,
    b: &Record,
) -> Option<FragFileFragment> {
    let ai = FragReadInfo::from(a);
    let bi = FragReadInfo::from(b);
    collect_fragment_with_frag_file_info(&ai, &bi)
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
        start: forward.pos,
        end: reverse.end,
        min_mapq: forward.mapq.min(reverse.mapq),
        read1_strand: if forward.is_read_1 {
            forward.strand
        } else {
            reverse.strand
        },
    })
}

/// Build a `FragFileFragment` from a single read (single-end input).
pub fn collect_fragment_with_frag_file_info_from_single_read(
    read: &FragReadInfo,
) -> Option<FragFileFragment> {
    if read.end <= read.pos {
        return None;
    }

    Some(FragFileFragment {
        tid: read.tid,
        start: read.pos,
        end: read.end,
        min_mapq: read.mapq,
        // No read 1/2 so just the strand of the single read
        read1_strand: read.strand,
    })
}
