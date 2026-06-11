use anyhow::{Context, Result};
use rust_htslib::bam::Record;

use crate::shared::fragment::read_order_fragment::{
    FragmentWithReadOrder, ReadOrderReadInfo, collect_fragment_with_read_order,
    collect_fragment_with_read_order_from_single_read,
};

use super::{InputItem, Pairer, PairingAdapter};

/// Pairing policy for minimal read-order fragments.
///
/// This is used for count-only passes where keeping BAM records would be unnecessary work. The
/// pairing rules delegate to `collect_fragment_with_read_order` so this path accepts and rejects
/// the same read pairs as the full-record fragment iterator used by output-producing commands.
pub(crate) struct WithReadOrderPairer;

impl Pairer for WithReadOrderPairer {
    type Read = ReadOrderReadInfo;
    type Output = FragmentWithReadOrder;

    /// Pair two reads into a minimal fragment span.
    ///
    /// Invalid or ambiguously marked pairs return `None`, exactly as the full-record pairing path
    /// does for output fragments.
    fn pair(&self, first: &Self::Read, second: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_read_order(first, second)
    }
}

/// Create a BAM fragment iterator that keeps only span and read-order checks.
///
/// The iterator applies the same read inclusion hook and fragment filter shape as
/// `fragments_with_records_from_bam`, but maps each BAM record into a lightweight read-order
/// record before pairing. Allelic-fragments uses this path for the control-count pass, so changes
/// here should be checked against the full-record iterator before changing which pairs are accepted.
pub(crate) fn fragments_with_read_order_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    fragment_filter: PF,
    unpaired: bool,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem>>,
    WithReadOrderPairer,
    ReadOrderReadInfo,
    FragmentWithReadOrder,
>
where
    RIter: Iterator<Item = Result<Record>>,
    PF: Fn(&FragmentWithReadOrder) -> bool + Send + Sync + 'static,
{
    let mapped = records.map(|result| {
        result
            .context("reading BAM record")
            .map(InputItem::BamRecord)
    });

    let mut adapter = PairingAdapter::new(
        mapped,
        if unpaired {
            None::<WithReadOrderPairer>
        } else {
            Some(WithReadOrderPairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, |record| {
        ReadOrderReadInfo::try_from(record).map_err(anyhow::Error::from)
    })
    .with_fragment_filter(fragment_filter);

    if unpaired {
        adapter = adapter
            .with_bam_single_fragment_from_read(collect_fragment_with_read_order_from_single_read);
    }

    adapter
}
