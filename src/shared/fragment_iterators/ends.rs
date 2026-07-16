use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::{
    commands::ends::config_structs::{BaseQualityFilter, ClipStrategy, KmerSource},
    shared::{
        fragment::ends_fragment::{
            EndReadInfo, FragmentWithEnds, collect_fragment_with_ends,
            collect_fragment_with_ends_from_single_read,
        },
        indel_mode::IndelMotifFilterPolicy,
    },
};

use super::{InputItem, Pairer, PairingAdapter};

#[derive(Clone)]
pub(crate) struct WithEndsPairer {
    pub(crate) clip_strategy: ClipStrategy,
    pub(crate) source_inside: KmerSource,
    pub(crate) indel_filter: IndelMotifFilterPolicy,
    pub(crate) k_inside: usize,
    pub(crate) max_soft_clips: u32,
    pub(crate) bq_filters: Vec<BaseQualityFilter>,
}

impl Pairer for WithEndsPairer {
    type Read = EndReadInfo;
    type Output = FragmentWithEnds;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_ends(
            a,
            b,
            self.clip_strategy,
            self.source_inside,
            self.indel_filter,
            self.k_inside,
            self.max_soft_clips,
            &self.bq_filters,
        )
    }
}

pub(crate) fn fragments_with_ends_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    clip_strategy: ClipStrategy,
    source_inside: KmerSource,
    indel_filter: IndelMotifFilterPolicy,
    k_inside: usize,
    max_soft_clips: u32,
    bq_filters: &[BaseQualityFilter],
    gc_tag: Option<&[u8]>,
    fragment_filter: PF,
    unpaired: bool,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem>>,
    WithEndsPairer,
    EndReadInfo,
    FragmentWithEnds,
>
where
    RIter: Iterator<Item = Result<Record>>,
    PF: Fn(&FragmentWithEnds) -> bool + Send + Sync + 'static,
{
    let pairer = WithEndsPairer {
        clip_strategy,
        source_inside,
        indel_filter,
        k_inside,
        max_soft_clips,
        bq_filters: bq_filters.to_vec(),
    };
    let gc_tag_bytes = gc_tag.map(|tag| tag.to_vec());
    let load_base_qualities = !bq_filters.is_empty();
    // Map BAM records -> InputItem::BamRecord, converting read errors to anyhow with context
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<WithEndsPairer>
        } else {
            Some(pairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, move |rec| {
        EndReadInfo::from_record_with_gc_tag(
            rec,
            gc_tag_bytes.as_deref(),
            clip_strategy,
            k_inside,
            load_base_qualities,
        )
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        let bq_filters = bq_filters.to_vec();
        adapter = adapter.with_bam_single_fragment_from_read(move |read| {
            collect_fragment_with_ends_from_single_read(
                read,
                clip_strategy,
                source_inside,
                indel_filter,
                k_inside,
                max_soft_clips,
                &bq_filters,
            )
        });
    }

    adapter
}
