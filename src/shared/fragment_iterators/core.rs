// Shared iterator adaptor for BAM pairing and pre-built fragments.
//
// ------------------------ Usage examples ------------------------
/*
Example 1 (BAM, local counters — fastest per-thread):

let include_read = |rec: &Record| true;
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
    false,             // unpaired
).with_local_counters();

for frag in it.by_ref() {
    let frag = frag?;
    // use frag
}
let stats = it.counters_snapshot();
eprintln!("local: {:?}", stats);

Example 2 (BAM, shared counters across multiple iterators/threads):

let shared = SharedCounters::new(); // Cloneable handle
let mut it1 = fragments_from_bam(
    r1.records().map(|r| r.map_err(anyhow::Error::from)),
    include_read_basic,
    None,
    post_basic,
    false,
).with_shared_counters(shared.clone());

let mut it2 = fragments_from_bam(
    r2.records().map(|r| r.map_err(anyhow::Error::from)),
    include_read_basic,
    None,
    post_basic,
    false,
).with_shared_counters(shared.clone());

for _ in it1 {}
for _ in it2 {}

let stats = shared.snapshot();
eprintln!("shared: {:?}", stats);

Example 3 (Ready-made fragments, still want yielded count):

let mut it = fragments_from_iter(frag_iter_anyhow, |f| {
    let len = f.len();
    len >= min && len <= max
}).with_local_counters();

for _ in it.by_ref() {}
let stats = it.counters_snapshot();
*/

use anyhow::{Result, anyhow};
use fxhash::FxHashMap;
use rust_htslib::bam::Record;

use crate::shared::{
    fragment::{minimal_fragment::MinimalReadInfo, segment_fragment::SegmentedReadInfo},
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
    last_bam_coord: Option<(i32, i64)>,
    fragment_filter: Option<Box<dyn Fn(&F) -> bool + Send + Sync>>,
    counters: Box<dyn FragmentCounters + Send>,
    bam_include_read: Option<Box<dyn Fn(&Record) -> bool + Send + Sync>>,
    bam_map_read: Option<Box<dyn Fn(&Record) -> Result<R> + Send + Sync>>,
    // Optional converter used only when `pairer` is None (unpaired --reads-are-fragments mode).
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
            last_bam_coord: None,
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
        map_read: impl Fn(&Record) -> Result<R> + Send + Sync + 'static,
    ) -> Self {
        self.bam_include_read = Some(Box::new(include_read));
        self.bam_map_read = Some(Box::new(map_read));
        self
    }

    /// Unpaired: set the mapped-read -> fragment converter (used only when `pairer` is `None`).
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

#[inline]
fn ensure_nondecreasing_bam_coordinates(
    last_bam_coord: &mut Option<(i32, i64)>,
    rec: &Record,
) -> Result<()> {
    let current = (rec.tid(), rec.pos());
    if let Some(previous) = *last_bam_coord {
        if current.0 != previous.0 {
            return Err(anyhow!(
                "BAM reader yielded records from multiple tids inside single-chromosome stream: observed tid {} after tid {}",
                current.0,
                previous.0
            ));
        }
        if current.1 < previous.1 {
            return Err(anyhow!(
                "BAM records must be coordinate-sorted with nondecreasing read.pos within single-chromosome stream, but observed pos {} after {} on tid {}",
                current.1,
                previous.1,
                current.0
            ));
        }
    }
    *last_bam_coord = Some(current);
    Ok(())
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
                    if let Some(accept_fragment) = &self.fragment_filter
                        && !accept_fragment(&f)
                    {
                        continue;
                    }
                    self.counters.inc_yielded_fragments();
                    return Some(Ok(f));
                }
                Ok(InputItem::BamRecord(rec)) => {
                    if let Err(error) =
                        ensure_nondecreasing_bam_coordinates(&mut self.last_bam_coord, &rec)
                    {
                        return Some(Err(error));
                    }
                    // Count every incoming BAM record
                    self.counters.inc_incoming_reads();
                    // Apply include_read if configured
                    if let Some(pred) = &self.bam_include_read
                        && !pred(&rec)
                    {
                        continue;
                    }
                    // Accepted read by initial flag / mapq policy
                    self.counters.inc_accepted_reads(rec.is_reverse());
                    let Some(map_read) = &self.bam_map_read else {
                        return Some(Err(anyhow!("BAM record seen but no mapper configured")));
                    };
                    let mapped = match map_read(&rec) {
                        Ok(mapped_read) => mapped_read,
                        Err(error) => {
                            return Some(Err(error.context("mapping BAM record")));
                        }
                    };

                    // Unpaired path when no pairer is present
                    if self.pairer.is_none() {
                        let Some(map_frag) = &self.bam_single_fragment_from_read else {
                            return Some(Err(anyhow!(
                                "no pairer and unpaired fragment mapper not configured"
                            )));
                        };
                        let frag_opt = map_frag(&mapped);
                        if let Some(frag) = frag_opt {
                            self.counters.inc_produced_fragments();
                            if let Some(accept_fragment) = &self.fragment_filter
                                && !accept_fragment(&frag)
                            {
                                continue;
                            }
                            self.counters.inc_yielded_fragments();
                            return Some(Ok(frag));
                        }
                        continue;
                    }

                    // Paired-end path: stash by qname and emit when both mates are available
                    let qname = rec.qname().to_vec();
                    let read = mapped;
                    if let Some(mate) = self.stash.remove(&qname) {
                        let Some(pairer) = self.pairer.as_ref() else {
                            return Some(Err(anyhow!("pairer required for BAM reads")));
                        };
                        if let Some(frag) = pairer.pair(&read, &mate) {
                            self.counters.inc_produced_fragments();
                            if let Some(accept_fragment) = &self.fragment_filter
                                && !accept_fragment(&frag)
                            {
                                continue;
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

#[cfg(test)]
mod tests {
    include!("core_tests.rs");
}
