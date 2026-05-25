use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::shared::{
    fragment::indel_counting_fragment::{
        FragmentWithIndelCounts, IndelReadInfo, collect_fragment_with_indel_counts,
        collect_fragment_with_indel_counts_from_single_read,
    },
    indel_mode::IndelMode,
};

use super::{InputItem, Pairer, PairingAdapter};

/* WithIndelCounts pairing */

pub(crate) struct WithIndelCountsPairer {
    pub(crate) indel_mode: IndelMode,
}

impl Pairer for WithIndelCountsPairer {
    type Read = IndelReadInfo;
    type Output = FragmentWithIndelCounts;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_indel_counts(
            a,
            b,
            matches!(self.indel_mode, IndelMode::Skip),
            matches!(self.indel_mode, IndelMode::Adjust),
        )
    }
}

pub(crate) type IndelCountsIter<'a> = PairingAdapter<
    Box<dyn Iterator<Item = Result<InputItem>> + 'a>,
    WithIndelCountsPairer,
    IndelReadInfo,
    FragmentWithIndelCounts,
>;

pub(crate) fn fragments_with_indel_counts_from_bam<'a, RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    indel_mode: IndelMode,
    inspect_cigar: bool,
    fragment_filter: PF,
    unpaired: bool,
) -> IndelCountsIter<'a>
where
    RIter: Iterator<Item = Result<Record>> + 'a,
    PF: Fn(&FragmentWithIndelCounts) -> bool + Send + Sync + 'static,
{
    let pairer = WithIndelCountsPairer { indel_mode };
    // Map BAM records -> InputItem::BamRecord, converting read errors to anyhow with context
    let mapped: Box<dyn Iterator<Item = Result<InputItem>> + 'a> =
        Box::new(records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord)));

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<WithIndelCountsPairer>
        } else {
            Some(pairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, move |rec| {
        IndelReadInfo::from_record(rec, inspect_cigar).map_err(anyhow::Error::from)
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        let skip_indels = matches!(indel_mode, IndelMode::Skip);
        let count_indels = matches!(indel_mode, IndelMode::Adjust);
        adapter = adapter.with_bam_single_fragment_from_read(move |read| {
            collect_fragment_with_indel_counts_from_single_read(read, skip_indels, count_indels)
        });
    }

    adapter
}
