use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::Record;

/// Compact per-read info with extracted indel events.
#[derive(Debug, Clone)]
pub struct WithRecordReadInfo {
    pub tid: i32,
    pub pos: u32, // Leftmost aligned reference pos
    pub end: u32, // Exclusive rightmost aligned reference end
    pub is_reverse: bool,
    pub mapq: u8,
    pub record: Record,
}

impl From<&Record> for WithRecordReadInfo {
    #[inline]
    fn from(r: &Record) -> Self {
        WithRecordReadInfo {
            tid: r.tid(),
            pos: r.pos() as u32,
            end: r.reference_end() as u32,
            is_reverse: r.is_reverse(),
            mapq: r.mapq(),
            record: r.clone(),
        }
    }
}

impl PairOrientable for WithRecordReadInfo {
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
pub struct WithRecordsFragment {
    pub tid: i32,
    pub start: u32, // forward.pos
    pub end: u32,   // reverse.end (end-exclusive)
    pub min_mapq: u8,
    pub single_record: Option<Record>,
    pub forward_record: Option<Record>,
    pub reverse_record: Option<Record>,
}

impl WithRecordsFragment {
    /// Reference-span fragment length (end - start).
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start
    }
}

/// Build a `WithRecordsFragment` from two `Record`s.
///
/// NOTE: Consumes the records.
#[inline]
pub fn collect_fragment_with_records_from_records(
    a: &Record,
    b: &Record,
) -> Option<WithRecordsFragment> {
    let ai = WithRecordReadInfo::from(a);
    let bi = WithRecordReadInfo::from(b);
    collect_fragment_with_records(&ai, &bi)
}

/// Build a `WithRecordsFragment` from two reads.
///
/// NOTE: Consumes the records.
pub fn collect_fragment_with_records(
    a: &WithRecordReadInfo,
    b: &WithRecordReadInfo,
) -> Option<WithRecordsFragment> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }

    Some(WithRecordsFragment {
        tid: forward.tid,
        start: forward.pos,
        end: reverse.end,
        min_mapq: forward.mapq.min(reverse.mapq),
        single_record: None,
        // TODO: Avoid cloning. Would like to keep reusing oriented_pair_from_read_info but perhaps an owned version of it is needed?
        forward_record: Some(forward.record.clone()),
        reverse_record: Some(reverse.record.clone()),
    })
}

/// Build a `WithRecordsFragment` from a single read (unpaired input).
pub fn collect_fragment_with_records_from_single_read(
    read: &WithRecordReadInfo,
) -> Option<WithRecordsFragment> {
    if read.end <= read.pos {
        return None;
    }

    Some(WithRecordsFragment {
        tid: read.tid,
        start: read.pos,
        end: read.end,
        min_mapq: read.mapq,
        single_record: Some(read.record.clone()),
        forward_record: None,
        reverse_record: None,
    })
}
