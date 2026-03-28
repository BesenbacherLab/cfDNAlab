use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::interval::Interval;
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::Record;

/// Compact per-read info with extracted indel events.
#[derive(Debug, Clone)]
pub struct WithRecordReadInfo {
    pub tid: i32,
    pub interval: Interval<u32>, // Aligned reference span [start: pos(), end: reference_end())
    pub is_reverse: bool,
    pub mapq: u8,
    pub record: Record,
}

impl TryFrom<&Record> for WithRecordReadInfo {
    type Error = crate::Error;

    #[inline]
    fn try_from(r: &Record) -> Result<Self> {
        Ok(WithRecordReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            mapq: r.mapq(),
            record: r.clone(),
        })
    }
}

impl WithRecordReadInfo {
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
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
        self.start()
    }
}

/// Fragment with mapq and read1 strand.
#[derive(Debug, Clone)]
pub struct WithRecordsFragment {
    pub tid: i32,
    pub interval: Interval<u32>, // forward.pos .. reverse.end
    pub min_mapq: u8,
    pub single_record: Option<Record>,
    pub forward_record: Option<Record>,
    pub reverse_record: Option<Record>,
}

impl WithRecordsFragment {
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

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
        interval: Interval::new(forward.start(), reverse.end()).ok()?,
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
    Some(WithRecordsFragment {
        tid: read.tid,
        interval: read.interval,
        min_mapq: read.mapq,
        single_record: Some(read.record.clone()),
        forward_record: None,
        reverse_record: None,
    })
}
