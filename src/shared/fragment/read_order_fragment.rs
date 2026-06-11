use crate::Result;
use crate::shared::fragment::minimal_fragment::{
    PairOrientable, is_inwards_oriented, oriented_pair_from_read_info,
};
use crate::shared::interval::Interval;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Record;

/// Minimal per-read info for the allelic-fragments control-counting pass.
///
/// This type keeps only the fields needed to reconstruct the project fragment span and paired-read
/// validity. It deliberately does not clone the BAM record, because pass 1 only counts case/control
/// candidates and never inspects bases.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadOrderReadInfo {
    pub(crate) tid: i32,
    pub(crate) interval: Interval<u32>,
    pub(crate) is_reverse: bool,
    pub(crate) is_read_1: bool,
}

impl TryFrom<&Record> for ReadOrderReadInfo {
    type Error = crate::Error;

    /// Build minimal read info from a BAM record.
    ///
    /// The interval stores `pos` to `reference_end`, matching the full-record fragment path. Any
    /// invalid interval is surfaced here so later pairing code only handles well-formed read spans.
    #[inline]
    fn try_from(record: &Record) -> Result<Self> {
        Ok(Self {
            tid: record.tid(),
            interval: Interval::new(record.pos() as u32, record.reference_end() as u32)?,
            is_reverse: record.is_reverse(),
            is_read_1: record.is_first_in_template(),
        })
    }
}

impl ReadOrderReadInfo {
    /// Aligned reference start.
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Aligned reference end.
    #[inline]
    pub(crate) fn end(&self) -> u32 {
        self.interval.end()
    }
}

impl PairOrientable for ReadOrderReadInfo {
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

/// Minimal fragment span used to count controls without cloning BAM records.
///
/// The span uses the same project semantics as `WithRecordsFragment`: paired fragments run from
/// `forward.pos` to `reverse.reference_end`, while unpaired reads use the aligned read interval.
#[derive(Debug, Clone, Copy)]
pub(crate) struct FragmentWithReadOrder {
    pub(crate) interval: Interval<u32>,
}

impl FragmentWithReadOrder {
    /// Fragment start.
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Fragment end.
    #[inline]
    pub(crate) fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Reference-span fragment length.
    #[inline]
    pub(crate) fn len(&self) -> u32 {
        self.interval.len()
    }
}

/// Build a paired fragment span without keeping the BAM records.
///
/// This must accept and reject the same read pairs as `collect_fragment_with_records`: both reads
/// must be on the same contig, face inward, and have exactly one read1 mate. Allelic-fragments uses
/// this path for counting controls and the record-retaining path for writing output, so a mismatch
/// would make the two passes disagree about which fragments exist.
pub(crate) fn collect_fragment_with_read_order(
    first: &ReadOrderReadInfo,
    second: &ReadOrderReadInfo,
) -> Option<FragmentWithReadOrder> {
    let (forward, reverse) = oriented_pair_from_read_info(first, second)?;
    if !is_inwards_oriented(forward, reverse) {
        return None;
    }
    if forward.is_read_1 == reverse.is_read_1 {
        return None;
    }

    Some(FragmentWithReadOrder {
        interval: Interval::new(forward.start(), reverse.end()).ok()?,
    })
}

/// Build a minimal fragment from one read in `--reads-are-fragments` mode.
///
/// The read-order count pass uses this for unpaired input so the counted fragment span is identical
/// to the full-record output pass. Read inclusion filtering has already removed reads that should
/// not be considered fragments.
pub(crate) fn collect_fragment_with_read_order_from_single_read(
    read: &ReadOrderReadInfo,
) -> Option<FragmentWithReadOrder> {
    Some(FragmentWithReadOrder {
        interval: read.interval,
    })
}
