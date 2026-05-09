use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::shared::fragment::segment_fragment::{
    FragmentWithSegments, SegmentedReadInfo, collect_fragment_with_segments,
    collect_fragment_with_segments_from_single_read,
};

use super::{InputItem, Pairer, PairingAdapter};

pub struct WithSegmentsPairer {
    pub trigger_min_gap_bp: u32,
    pub include_inter_mate_gap: bool,
}

impl Pairer for WithSegmentsPairer {
    type Read = SegmentedReadInfo;
    type Output = FragmentWithSegments;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_segments(a, b, self.trigger_min_gap_bp, self.include_inter_mate_gap)
    }
}

pub fn fragments_with_segments_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    trigger_min_gap_bp: u32,
    include_inter_mate_gap: bool,
    gc_tag: Option<&[u8]>,
    fragment_filter: PF,
    unpaired: bool,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem<FragmentWithSegments>>>,
    WithSegmentsPairer,
    SegmentedReadInfo,
    FragmentWithSegments,
>
where
    RIter: Iterator<Item = Result<Record>>,
    PF: Fn(&FragmentWithSegments) -> bool + Send + Sync + 'static,
{
    let pairer = WithSegmentsPairer {
        trigger_min_gap_bp,
        include_inter_mate_gap,
    };
    let gc_tag_bytes = gc_tag.map(|tag| tag.to_vec());
    // Map BAM records -> InputItem::BamRecord, converting read errors to anyhow with context
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<WithSegmentsPairer>
        } else {
            Some(pairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, move |rec| {
        SegmentedReadInfo::from_record_with_gc_tag(rec, gc_tag_bytes.as_deref())
            .map_err(anyhow::Error::from)
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        adapter = adapter.with_bam_single_fragment_from_read(move |read| {
            collect_fragment_with_segments_from_single_read(read, trigger_min_gap_bp)
        });
    }

    adapter
}

pub fn fragments_with_segments_from_iter<I, PF>(
    frags: I,
    fragment_filter: PF,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem<FragmentWithSegments>>>,
    WithSegmentsPairer,
    SegmentedReadInfo,
    FragmentWithSegments,
>
where
    I: Iterator<Item = Result<FragmentWithSegments>>,
    PF: Fn(&FragmentWithSegments) -> bool + Send + Sync + 'static,
{
    let mapped = frags.map(|res| res.map(InputItem::Fragment));
    PairingAdapter::new(mapped, None::<WithSegmentsPairer>).with_fragment_filter(fragment_filter)
}
