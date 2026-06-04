use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::shared::fragment::frag_file_fragment::{
    FragFileFragment, FragReadInfo, collect_fragment_with_frag_file_info,
    collect_fragment_with_frag_file_info_from_single_read,
};

use super::{InputItem, Pairer, PairingAdapter};

/* For frag files pairing */

pub(crate) struct WithFragInfoPairer;

impl Pairer for WithFragInfoPairer {
    type Read = FragReadInfo;
    type Output = FragFileFragment;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_frag_file_info(a, b)
    }
}

pub(crate) fn fragments_with_frag_file_info_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    fragment_filter: PF,
    unpaired: bool,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem>>,
    WithFragInfoPairer,
    FragReadInfo,
    FragFileFragment,
>
where
    RIter: Iterator<Item = Result<Record>>,
    PF: Fn(&FragFileFragment) -> bool + Send + Sync + 'static,
{
    // Map BAM records -> InputItem::BamRecord, converting read errors to anyhow with context
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<WithFragInfoPairer>
        } else {
            Some(WithFragInfoPairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, |rec| {
        FragReadInfo::try_from(rec).map_err(anyhow::Error::from)
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        adapter = adapter.with_bam_single_fragment_from_read(|read| {
            collect_fragment_with_frag_file_info_from_single_read(read)
        });
    }

    adapter
}
