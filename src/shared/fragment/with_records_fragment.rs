use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::interval::Interval;
use rust_htslib::bam::ext::BamRecordExtensions; // reference_end()
use rust_htslib::bam::record::Record;

/// Per-read info used when fragment output needs the original BAM record.
///
/// This keeps the orientation fields needed for pairing plus the cloned `Record`. The
/// paired-fragment validity rules in this module must stay aligned with
/// `collect_fragment_with_read_order`, because allelic-fragments uses this representation for row
/// output and the minimal read-order representation for the control-count pass.
#[derive(Debug, Clone)]
pub(crate) struct WithRecordReadInfo {
    pub(crate) tid: i32,
    pub(crate) interval: Interval<u32>, // Aligned reference span [start: pos(), end: reference_end())
    pub(crate) is_reverse: bool,
    pub(crate) record: Record,
}

impl TryFrom<&Record> for WithRecordReadInfo {
    type Error = crate::Error;

    /// Build full-record read info from a BAM record.
    ///
    /// The interval stores `pos` to `reference_end`.
    #[inline]
    fn try_from(r: &Record) -> Result<Self> {
        Ok(WithRecordReadInfo {
            tid: r.tid(),
            interval: Interval::new(r.pos() as u32, r.reference_end() as u32)?,
            is_reverse: r.is_reverse(),
            record: r.clone(),
        })
    }
}

impl WithRecordReadInfo {
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub(crate) fn end(&self) -> u32 {
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

/// Fragment with the BAM record or records used to build it.
#[derive(Debug, Clone)]
pub(crate) struct WithRecordsFragment {
    pub(crate) interval: Interval<u32>, // forward.pos .. reverse.end
    pub(crate) single_record: Option<Record>,
    pub(crate) forward_record: Option<Record>,
    pub(crate) reverse_record: Option<Record>,
}

impl WithRecordsFragment {
    #[inline]
    #[allow(dead_code)]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    #[allow(dead_code)]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Reference-span fragment length (end - start).
    #[inline]
    pub(crate) fn len(&self) -> u32 {
        self.interval.len()
    }
}

/// Build a `WithRecordsFragment` from two reads.
///
/// The reads must be on the same contig, inward oriented, and exactly one read must be marked
/// read1. The read1/read2 check rejects duplicate-mate or ambiguous pairs before row construction.
///
/// NOTE: Consumes the records.
pub(crate) fn collect_fragment_with_records(
    a: &WithRecordReadInfo,
    b: &WithRecordReadInfo,
) -> Option<WithRecordsFragment> {
    let (forward, reverse) = oriented_pair_from_read_info(a, b)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }
    if forward.record.is_first_in_template() == reverse.record.is_first_in_template() {
        return None;
    }

    Some(WithRecordsFragment {
        interval: Interval::new(forward.start(), reverse.end()).ok()?,
        single_record: None,
        // TODO: Avoid cloning. Would like to keep reusing oriented_pair_from_read_info but perhaps an owned version of it is needed?
        forward_record: Some(forward.record.clone()),
        reverse_record: Some(reverse.record.clone()),
    })
}

/// Build a `WithRecordsFragment` from a single read in `--reads-are-fragments` mode.
///
/// Read filtering has already decided whether the record is acceptable. No mate-orientation checks
/// are applied because the read itself is the fragment.
pub(crate) fn collect_fragment_with_records_from_single_read(
    read: &WithRecordReadInfo,
) -> Option<WithRecordsFragment> {
    Some(WithRecordsFragment {
        interval: read.interval,
        single_record: Some(read.record.clone()),
        forward_record: None,
        reverse_record: None,
    })
}
