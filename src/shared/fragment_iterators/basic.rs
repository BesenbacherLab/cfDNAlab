use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::shared::fragment::minimal_fragment::{
    Fragment, MinimalReadInfo, collect_fragment, collect_fragment_from_single_read,
};

use super::{InputItem, Pairer, PairingAdapter};

pub(crate) struct BasicPairer;

impl Pairer for BasicPairer {
    type Read = MinimalReadInfo;
    type Output = Fragment;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment(a, b)
    }
}

pub(crate) type BasicFragmentIter<'a> = PairingAdapter<
    Box<dyn Iterator<Item = Result<InputItem>> + 'a>,
    BasicPairer,
    MinimalReadInfo,
    Fragment,
>;

pub(crate) fn fragments_from_bam<'a, RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    gc_tag: Option<&[u8]>,
    fragment_filter: PF,
    unpaired: bool,
) -> BasicFragmentIter<'a>
where
    RIter: Iterator<Item = Result<Record>> + 'a,
    PF: Fn(&Fragment) -> bool + Send + Sync + 'static,
{
    let gc_tag_bytes = gc_tag.map(|tag| tag.to_vec());
    let mapped: Box<dyn Iterator<Item = Result<InputItem>> + 'a> =
        Box::new(records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord)));

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<BasicPairer>
        } else {
            Some(BasicPairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, move |rec| {
        MinimalReadInfo::from_record_with_gc_tag(rec, gc_tag_bytes.as_deref())
            .map_err(anyhow::Error::from)
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        adapter = adapter.with_bam_single_fragment_from_read(collect_fragment_from_single_read);
    }

    adapter
}
