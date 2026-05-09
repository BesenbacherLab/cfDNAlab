use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::shared::fragment::with_records_fragment::{
    WithRecordReadInfo, WithRecordsFragment, collect_fragment_with_records,
    collect_fragment_with_records_from_single_read,
};

use super::{InputItem, Pairer, PairingAdapter};

pub struct WithRecordReadInfoPairer;

impl Pairer for WithRecordReadInfoPairer {
    type Read = WithRecordReadInfo;
    type Output = WithRecordsFragment;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_records(a, b)
    }
}

pub fn fragments_with_records_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    fragment_filter: PF,
    unpaired: bool,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem<WithRecordsFragment>>>,
    WithRecordReadInfoPairer,
    WithRecordReadInfo,
    WithRecordsFragment,
>
where
    RIter: Iterator<Item = Result<Record>>,
    PF: Fn(&WithRecordsFragment) -> bool + Send + Sync + 'static,
{
    // Map BAM records -> InputItem::BamRecord, converting read errors to anyhow with context
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<WithRecordReadInfoPairer>
        } else {
            Some(WithRecordReadInfoPairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, |rec| {
        WithRecordReadInfo::try_from(rec).map_err(anyhow::Error::from)
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        adapter = adapter.with_bam_single_fragment_from_read(|read| {
            collect_fragment_with_records_from_single_read(read)
        });
    }

    adapter
}
