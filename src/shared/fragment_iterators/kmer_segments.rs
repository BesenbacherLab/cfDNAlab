use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::shared::{
    fragment::segment_kmer_fragment::{
        FragmentWithKmerSegments, KmerSegmentedReadInfo, collect_fragment_with_kmer_segments,
        collect_fragment_with_kmer_segments_from_single_read,
    },
    indel_mode::IndelMode,
};

use super::{InputItem, Pairer, PairingAdapter};

/* Kmer segments pairing */

#[derive(Clone, Copy)]
pub(crate) struct KmerSegmentsPairer {
    pub(crate) indel_mode: IndelMode,
    pub(crate) include_inter_mate_gap: bool,
    pub(crate) end_offset: u32,
}

impl Pairer for KmerSegmentsPairer {
    type Read = KmerSegmentedReadInfo;
    type Output = FragmentWithKmerSegments;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_kmer_segments(
            a,
            b,
            self.indel_mode,
            self.include_inter_mate_gap,
            self.end_offset,
        )
    }
}

pub(crate) fn fragments_with_kmer_segments_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    indel_mode: IndelMode,
    include_inter_mate_gap: bool,
    end_offset: u32,
    gc_tag: Option<&[u8]>,
    fragment_filter: PF,
    unpaired: bool,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem>>,
    KmerSegmentsPairer,
    KmerSegmentedReadInfo,
    FragmentWithKmerSegments,
>
where
    RIter: Iterator<Item = Result<Record>>,
    PF: Fn(&FragmentWithKmerSegments) -> bool + Send + Sync + 'static,
{
    let pairer = KmerSegmentsPairer {
        indel_mode,
        include_inter_mate_gap,
        end_offset,
    };
    let gc_tag_bytes = gc_tag.map(|tag| tag.to_vec());
    // Map BAM records -> InputItem::BamRecord, converting read errors to anyhow with context
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));
    let capture_segments = matches!(indel_mode, IndelMode::Adjust);

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<KmerSegmentsPairer>
        } else {
            Some(pairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, move |rec| {
        KmerSegmentedReadInfo::from_record(rec, capture_segments, gc_tag_bytes.as_deref())
            .map_err(anyhow::Error::from)
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        adapter = adapter.with_bam_single_fragment_from_read(move |read| {
            collect_fragment_with_kmer_segments_from_single_read(read, indel_mode, end_offset)
        });
    }

    adapter
}
