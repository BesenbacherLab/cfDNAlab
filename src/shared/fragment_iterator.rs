// This module unifies BAM-based pairing and ready-made fragments under
// a single iterator interface, with pluggable pairing policies.

// ------------------------ Usage examples ------------------------
/*
Example 1 (BAM, local counters — fastest per-thread):

let include_read = |rec: &Record| /* your flag policy */ true;
let post = |f: &FragmentWithSegments| {
    let len = f.len();
    (len >= opt.fragment_lengths.min_fragment_length) &&
    (len <= opt.fragment_lengths.max_fragment_length)
};

let mut it = fragments_with_segments_from_bam(
    reader.records().map(|r| r.map_err(anyhow::Error::from)),
    include_read,
    1,                 // trigger_min_gap_bp
    !opt.ignore_gap,   // include_inter_mate_gap
    None,              // optional gc tag
    post,
).with_local_counters();

for frag in it.by_ref() {
    let frag = frag?;
    // use frag
}
let stats = it.counters_snapshot();
eprintln!("local: {:?}", stats);

Example 2 (BAM, **shared** counters across multiple iterators/threads):

let shared = SharedCounters::new(); // Cloneable handle
let mut it1 = fragments_basic_from_bam(
    r1.records().map(|r| r.map_err(anyhow::Error::from)),
    include_read_basic,
    post_basic,
).with_shared_counters(shared.clone());

let mut it2 = fragments_basic_from_bam(
    r2.records().map(|r| r.map_err(anyhow::Error::from)),
    include_read_basic,
    post_basic,
).with_shared_counters(shared.clone());

// Run iterators (possibly on different threads)...
for _ in it1 {}
for _ in it2 {}

// Snapshot totals from the shared handle:
let stats = shared.snapshot();
eprintln!("shared: {:?}", stats);

Example 3 (Frag file iterator, no pairing, still want yielded count):

let mut it = fragments_basic_from_iter(frag_iter_anyhow, |f| {
    let len = f.len();
    len >= min && len <= max
}).with_local_counters();

for _ in it.by_ref() {}
let stats = it.counters_snapshot();
*/

use anyhow::{Context, Result, anyhow};
use fxhash::FxHashMap;
use rust_htslib::bam::Record;

use crate::shared::{
    fragment::{
        frag_file_fragment::{
            FragFileFragment, FragReadInfo, collect_fragment_with_frag_file_info,
        },
        indel_counting_fragment::{
            FragmentWithIndelCounts, IndelReadInfo, collect_fragment_with_indel_counts,
            collect_fragment_with_indel_counts_from_single_read,
        },
        minimal_fragment::{
            Fragment, MinimalReadInfo, collect_fragment, collect_fragment_from_single_read,
        },
        segment_fragment::{
            FragmentWithSegments, SegmentedReadInfo, collect_fragment_with_segments,
        },
        segment_kmer_fragment::{
            FragmentWithKmerSegments, KmerSegmentedReadInfo, collect_fragment_with_kmer_segments,
        },
        with_records_fragment::{
            WithRecordReadInfo, WithRecordsFragment, collect_fragment_with_records,
        },
    },
    indel_mode::IndelMode,
    iterator_counter::{
        FragmentCounterSnapshot, FragmentCounters, LocalCounters, NoopCounters, SharedCounters,
    },
};

pub trait HasStrand {
    fn is_reverse(&self) -> bool;
}

/// Shared trait to expose counter snapshots on boxed iterators.
pub trait FragmentIterCounters {
    fn counters_snapshot(&self) -> FragmentCounterSnapshot;
}

impl HasStrand for SegmentedReadInfo {
    #[inline]
    fn is_reverse(&self) -> bool {
        self.is_reverse
    }
}

impl HasStrand for MinimalReadInfo {
    #[inline]
    fn is_reverse(&self) -> bool {
        self.is_reverse
    }
}

/// Normalized items flowing into the adaptor: either a read (paired later by qname),
/// or a ready-made fragment that passes through as-is.
pub enum InputItem<F> {
    BamRecord(Record),
    Fragment(F),
}

/// Policy for turning two reads into a fragment.
pub trait Pairer {
    type Read;
    type Output;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output>;
}

/// Iterator adaptor that consumes `InputItem`s, pairs reads by qname,
/// and yields fragments.
pub struct PairingAdapter<I, P, R, F> {
    inner: I,
    pairer: Option<P>,
    stash: FxHashMap<Vec<u8>, R>,
    fragment_filter: Option<Box<dyn Fn(&F) -> bool + Send + Sync>>,
    counters: Box<dyn FragmentCounters + Send>,
    bam_include_read: Option<Box<dyn Fn(&Record) -> bool + Send + Sync>>,
    bam_map_read: Option<Box<dyn Fn(&Record) -> R + Send + Sync>>,
    // Optional converter used only when `pairer` is None (single-end mode).
    bam_single_fragment_from_read: Option<Box<dyn Fn(&R) -> Option<F> + Send + Sync>>,
}

impl<I, P, R, F> PairingAdapter<I, P, R, F>
where
    I: Iterator<Item = Result<InputItem<F>>>,
    P: Pairer<Read = R, Output = F>,
{
    pub fn new(inner: I, pairer: Option<P>) -> Self {
        Self {
            inner,
            pairer,
            stash: FxHashMap::default(),
            fragment_filter: None,
            counters: Box::new(NoopCounters),
            bam_include_read: None,
            bam_map_read: None,
            bam_single_fragment_from_read: None,
        }
    }

    pub fn with_fragment_filter(mut self, f: impl Fn(&F) -> bool + Send + Sync + 'static) -> Self {
        self.fragment_filter = Some(Box::new(f));
        self
    }

    /// BAM-only: set the include_read predicate and the Record->R mapper.
    pub fn with_bam_filter_and_mapper(
        mut self,
        include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
        map_read: impl Fn(&Record) -> R + Send + Sync + 'static,
    ) -> Self {
        self.bam_include_read = Some(Box::new(include_read));
        self.bam_map_read = Some(Box::new(map_read));
        self
    }

    /// Single-end: set the mapped-read -> fragment converter (used only when `pairer` is `None`).
    pub fn with_bam_single_fragment_from_read(
        mut self,
        map_fragment: impl Fn(&R) -> Option<F> + Send + Sync + 'static,
    ) -> Self {
        self.bam_single_fragment_from_read = Some(Box::new(map_fragment));
        self
    }

    /// Overwrite counters with fast, thread-local counters.
    #[inline]
    pub fn with_local_counters(mut self) -> Self {
        self.counters = Box::new(LocalCounters::new());
        self
    }

    /// Overwrite counters with a shared, atomic implementation.
    /// Hold on to a clone of `shared` if you want to read totals externally.
    #[inline]
    pub fn with_shared_counters(mut self, shared: SharedCounters) -> Self {
        self.counters = Box::new(shared);
        self
    }

    /// Read counters at any time (e.g., after a `for` loop using `.by_ref()`).
    #[inline]
    pub fn counters_snapshot(&self) -> FragmentCounterSnapshot {
        self.counters.snapshot()
    }
}

// TODO: In tools like fcoverage where we use extra fetch halos, we might end up
// counting (stats counters) fragments that *fall just outside the tile cores*
// in multiple tiles! That means we cannot use the stats to say how many reads
// and fragments were actually present (almost but not completely)
// We should look into fixing this (although low priority)

impl<I, P, R, F> Iterator for PairingAdapter<I, P, R, F>
where
    I: Iterator<Item = Result<InputItem<F>>>,
    P: Pairer<Read = R, Output = F>,
{
    type Item = Result<F>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let next_in = self.inner.next()?;
            match next_in {
                Err(e) => return Some(Err(e)),
                Ok(InputItem::Fragment(f)) => {
                    // Fragment already assembled upstream
                    self.counters.inc_incoming_fragments();
                    if let Some(accept_fragment) = &self.fragment_filter {
                        if !accept_fragment(&f) {
                            continue;
                        }
                    }
                    self.counters.inc_yielded_fragments();
                    return Some(Ok(f));
                }
                Ok(InputItem::BamRecord(rec)) => {
                    // Count every incoming BAM record
                    self.counters.inc_incoming_reads();
                    // Apply include_read if configured
                    if let Some(pred) = &self.bam_include_read {
                        if !pred(&rec) {
                            continue;
                        }
                    }
                    // Accepted read (by initial flags)
                    self.counters.inc_accepted_reads(rec.is_reverse());
                    let Some(map_read) = &self.bam_map_read else {
                        return Some(Err(anyhow!("BAM record seen but no mapper configured")));
                    };
                    let mapped = map_read(&rec);

                    // *Single-end* path when no pairer is present
                    if self.pairer.is_none() {
                        let Some(map_frag) = &self.bam_single_fragment_from_read else {
                            return Some(Err(anyhow!("single-end fragment mapper not configured")));
                        };
                        let frag_opt = map_frag(&mapped);
                        if let Some(frag) = frag_opt {
                            self.counters.inc_produced_fragments();
                            if let Some(accept_fragment) = &self.fragment_filter {
                                if !accept_fragment(&frag) {
                                    continue;
                                }
                            }
                            self.counters.inc_yielded_fragments();
                            return Some(Ok(frag));
                        }
                        continue;
                    }

                    // *Paired-end* path: stash by qname and emit when both mates are available
                    let qname = rec.qname().to_vec();
                    let read = mapped;
                    if let Some(mate) = self.stash.remove(&qname) {
                        let Some(pairer) = self.pairer.as_ref() else {
                            return Some(Err(anyhow!("pairer required for BAM reads")));
                        };
                        if let Some(frag) = pairer.pair(&read, &mate) {
                            self.counters.inc_produced_fragments();
                            if let Some(accept_fragment) = &self.fragment_filter {
                                if !accept_fragment(&frag) {
                                    continue;
                                }
                            }
                            self.counters.inc_yielded_fragments();
                            return Some(Ok(frag));
                        } else {
                            continue;
                        }
                    } else {
                        self.stash.insert(qname, read);
                        continue;
                    }
                }
            }
        }
    }
}

impl<I, P, R, F> FragmentIterCounters for PairingAdapter<I, P, R, F>
where
    I: Iterator<Item = Result<InputItem<F>>>,
    P: Pairer<Read = R, Output = F>,
{
    #[inline]
    fn counters_snapshot(&self) -> FragmentCounterSnapshot {
        PairingAdapter::counters_snapshot(self)
    }
}

/* WithSegments pairing */

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

/// From BAM: pair reads into `FragmentWithSegments`.
pub fn fragments_with_segments_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    trigger_min_gap_bp: u32,
    include_inter_mate_gap: bool,
    gc_tag: Option<&[u8]>,
    fragment_filter: PF,
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

    // Map BAM records -> InputItem::Read, converting errors to anyhow with context.
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    PairingAdapter::new(mapped, Some(pairer))
        .with_bam_filter_and_mapper(include_read, move |rec| {
            SegmentedReadInfo::from_record_with_gc_tag(rec, gc_tag_bytes.as_deref())
        })
        .with_fragment_filter(fragment_filter)
}

/// From an iterator of ready-made `FragmentWithSegments` (e.g., BED-like source).
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

/* Kmer segments pairing */

#[derive(Clone, Copy)]
pub struct KmerSegmentsPairer {
    pub indel_mode: IndelMode,
    pub include_inter_mate_gap: bool,
    pub end_offset: u32,
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

pub fn fragments_with_kmer_segments_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    indel_mode: IndelMode,
    include_inter_mate_gap: bool,
    end_offset: u32,
    gc_tag: Option<&[u8]>,
    fragment_filter: PF,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem<FragmentWithKmerSegments>>>,
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

    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    let capture_segments = matches!(indel_mode, IndelMode::Adjust);

    PairingAdapter::new(mapped, Some(pairer))
        .with_bam_filter_and_mapper(include_read, move |rec| {
            KmerSegmentedReadInfo::from_record(rec, capture_segments, gc_tag_bytes.as_deref())
        })
        .with_fragment_filter(fragment_filter)
}

/* Basic fragment pairing */

pub struct BasicPairer;

impl Pairer for BasicPairer {
    type Read = MinimalReadInfo;
    type Output = Fragment;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment(a, b)
    }
}

// Shared alias to keep basic fragment iterators uniform across constructors
pub type BasicFragmentIter<'a> = PairingAdapter<
    Box<dyn Iterator<Item = Result<InputItem<Fragment>>> + 'a>,
    BasicPairer,
    MinimalReadInfo,
    Fragment,
>;

/// From BAM: optionally pair reads into `Fragment` or treat each as single-end.
pub fn fragments_from_bam<'a, RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    gc_tag: Option<&[u8]>,
    fragment_filter: PF,
    single_end: bool,
) -> BasicFragmentIter<'a>
where
    RIter: Iterator<Item = Result<Record>> + 'a,
    PF: Fn(&Fragment) -> bool + Send + Sync + 'static,
{
    let gc_tag_bytes = gc_tag.map(|tag| tag.to_vec());
    let mapped: Box<dyn Iterator<Item = Result<InputItem<Fragment>>> + 'a> =
        Box::new(records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord)));

    let mut adapter = PairingAdapter::new(
        mapped,
        if single_end {
            None::<BasicPairer>
        } else {
            Some(BasicPairer)
        },
    )
    .with_bam_filter_and_mapper(include_read, move |rec| {
        let mut info = MinimalReadInfo::from(rec);
        if let Some(tag) = gc_tag_bytes.as_deref() {
            info.gc_tag = crate::shared::gc_tag::read_gc_tag_from_record(rec, tag);
        }
        info
    })
    .with_fragment_filter(fragment_filter);

    if single_end {
        adapter = adapter
            .with_bam_single_fragment_from_read(|read| collect_fragment_from_single_read(read));
    }

    adapter
}

/// From an iterator of ready-made basic `Fragment`s.
pub fn fragments_from_iter<I, PF>(
    frags: I,
    fragment_filter: PF,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem<Fragment>>>,
    BasicPairer,
    MinimalReadInfo,
    Fragment,
>
where
    I: Iterator<Item = Result<Fragment>>,
    PF: Fn(&Fragment) -> bool + Send + Sync + 'static,
{
    let mapped = frags.map(|res| res.map(InputItem::Fragment));

    PairingAdapter::new(mapped, None::<BasicPairer>).with_fragment_filter(fragment_filter)
}

/* WithIndelCounts pairing */

pub struct WithIndelCountsPairer {
    pub indel_mode: IndelMode,
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

// Shared alias to keep fragment iterators uniform across constructors
pub type IndelCountsIter<'a> = PairingAdapter<
    Box<dyn Iterator<Item = Result<InputItem<FragmentWithIndelCounts>>> + 'a>,
    WithIndelCountsPairer,
    IndelReadInfo,
    FragmentWithIndelCounts,
>;

/// From BAM: pair reads into `FragmentWithIndelCounts`.
pub fn fragments_with_indel_counts_from_bam<'a, RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    indel_mode: IndelMode,
    fragment_filter: PF,
) -> IndelCountsIter<'a>
where
    RIter: Iterator<Item = Result<Record>> + 'a,
    PF: Fn(&FragmentWithIndelCounts) -> bool + Send + Sync + 'static,
{
    let pairer = WithIndelCountsPairer { indel_mode };

    // Map BAM records -> InputItem::Read, converting errors to anyhow with context.
    let mapped: Box<dyn Iterator<Item = Result<InputItem<FragmentWithIndelCounts>>> + 'a> =
        Box::new(records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord)));

    PairingAdapter::new(mapped, Some(pairer))
        .with_bam_filter_and_mapper(include_read, |rec| IndelReadInfo::from(rec))
        .with_fragment_filter(fragment_filter)
}

/// From an iterator of ready-made `FragmentWithIndelCounts` (e.g., BED-like source).
pub fn fragments_with_indel_counts_from_iter<'a, I, PF>(
    frags: I,
    fragment_filter: PF,
) -> IndelCountsIter<'a>
where
    I: Iterator<Item = Result<FragmentWithIndelCounts>> + 'a,
    PF: Fn(&FragmentWithIndelCounts) -> bool + Send + Sync + 'static,
{
    let mapped: Box<dyn Iterator<Item = Result<InputItem<FragmentWithIndelCounts>>> + 'a> =
        Box::new(frags.map(|res| res.map(InputItem::Fragment)));

    PairingAdapter::new(mapped, None::<WithIndelCountsPairer>).with_fragment_filter(fragment_filter)
}

/// From BAM: treat each read as a single-end fragment with optional indel adjustment.
pub fn fragments_with_indel_counts_from_single_end_bam<'a, RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    indel_mode: IndelMode,
    fragment_filter: PF,
) -> IndelCountsIter<'a>
where
    RIter: Iterator<Item = Result<Record>> + 'a,
    PF: Fn(&FragmentWithIndelCounts) -> bool + Send + Sync + 'static,
{
    let skip_indels = matches!(indel_mode, IndelMode::Skip);
    let count_indels = matches!(indel_mode, IndelMode::Adjust);

    // Map BAM records -> InputItem::BamRecord, converting errors to anyhow with context.
    let mapped: Box<dyn Iterator<Item = Result<InputItem<FragmentWithIndelCounts>>> + 'a> =
        Box::new(records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord)));

    PairingAdapter::new(mapped, None::<WithIndelCountsPairer>)
        .with_bam_filter_and_mapper(include_read, |rec| IndelReadInfo::from(rec))
        .with_bam_single_fragment_from_read(move |read| {
            collect_fragment_with_indel_counts_from_single_read(read, skip_indels, count_indels)
        })
        .with_fragment_filter(fragment_filter)
}

/* For frag files pairing */

pub struct WithFragInfoPairer;

impl Pairer for WithFragInfoPairer {
    type Read = FragReadInfo;
    type Output = FragFileFragment;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_frag_file_info(a, b)
    }
}

/// From BAM: pair reads into `FragFileFragment`.
pub fn fragments_with_frag_file_info_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    fragment_filter: PF,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem<FragFileFragment>>>,
    WithFragInfoPairer,
    FragReadInfo,
    FragFileFragment,
>
where
    RIter: Iterator<Item = Result<Record>>,
    PF: Fn(&FragFileFragment) -> bool + Send + Sync + 'static,
{
    let pairer = WithFragInfoPairer {};

    // Map BAM records -> InputItem::Read, converting errors to anyhow with context.
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    PairingAdapter::new(mapped, Some(pairer))
        .with_bam_filter_and_mapper(include_read, |rec| FragReadInfo::from(rec))
        .with_fragment_filter(fragment_filter)
}

/// From an iterator of ready-made `FragFileFragment` (e.g., BED-like source).
pub fn fragments_with_frag_file_info_from_iter<I, PF>(
    frags: I,
    fragment_filter: PF,
) -> PairingAdapter<
    impl Iterator<Item = Result<InputItem<FragFileFragment>>>,
    WithFragInfoPairer,
    FragReadInfo,
    FragFileFragment,
>
where
    I: Iterator<Item = Result<FragFileFragment>>,
    PF: Fn(&FragFileFragment) -> bool + Send + Sync + 'static,
{
    let mapped = frags.map(|res| res.map(InputItem::Fragment));

    PairingAdapter::new(mapped, None::<WithFragInfoPairer>).with_fragment_filter(fragment_filter)
}

/* For with-records read pairing */

pub struct WithRecordReadInfoPairer;

impl Pairer for WithRecordReadInfoPairer {
    type Read = WithRecordReadInfo;
    type Output = WithRecordsFragment;

    fn pair(&self, a: &Self::Read, b: &Self::Read) -> Option<Self::Output> {
        collect_fragment_with_records(a, b)
    }
}

/// From BAM: pair reads into `FragFileFragment`.
pub fn fragments_with_records_from_bam<RIter, PF>(
    records: RIter,
    include_read: impl Fn(&Record) -> bool + Send + Sync + 'static,
    fragment_filter: PF,
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
    let pairer = WithRecordReadInfoPairer {};

    // Map BAM records -> InputItem::Read, converting errors to anyhow with context.
    let mapped = records.map(|res| res.context("reading BAM record").map(InputItem::BamRecord));

    PairingAdapter::new(mapped, Some(pairer))
        .with_bam_filter_and_mapper(include_read, |rec| WithRecordReadInfo::from(rec))
        .with_fragment_filter(fragment_filter)
}

// /// From an iterator of ready-made `WithRecordsFragment` (e.g., BED-like source).
// pub fn fragments_with_records_from_iter<I, PF>(
//     frags: I,
//     fragment_filter: PF,
// ) -> PairingAdapter<
//     impl Iterator<Item = Result<InputItem<WithRecordsFragment>>>,
//     WithRecordReadInfoPairer,
//     WithRecordReadInfo,
//     WithRecordsFragment,
// >
// where
//     I: Iterator<Item = Result<WithRecordsFragment>>,
//     PF: Fn(&WithRecordsFragment) -> bool + Send + Sync + 'static,
// {
//     let mapped = frags.map(|res| res.map(InputItem::Fragment));

//     PairingAdapter::new(mapped, None::<WithRecordReadInfoPairer>)
//         .with_fragment_filter(fragment_filter)
// }
