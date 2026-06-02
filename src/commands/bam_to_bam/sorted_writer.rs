use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::sync::Arc;

use crate::shared::{
    constants::{
        COVERAGE_WEIGHT_AUX_TAG, FRAGMENT_COUNT_WEIGHT_AUX_TAG, FRAGMENT_LENGTH_AUX_TAG,
        GC_WEIGHT_AUX_TAG,
    },
    interval::Interval,
};
use anyhow::{Context, Result};
use rust_htslib::bam::{self, Record, record::Aux};

/// Per-record AUX tag data.
#[derive(Debug, Default)]
pub struct RecordTags {
    pub fragment_length: u32,
    pub coverage_weight: Option<f32>,
    pub fragment_count_weight: Option<f32>,
    pub gc_weight: Option<f32>,
}

/// Shared tag data for a pair of mate records.
pub type SharedTags = Arc<RecordTags>;

/// A buffered BAM record with its sort keys.
pub struct RecordEntry {
    pub interval: Interval<u32>,
    pub record: Record,
    pub tags: SharedTags,
}

impl RecordEntry {
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }
}

/// Keeps a bounded in-memory buffer so fragments are written in strict coordinate order
/// from an almost-sorted stream whose disorder is limited by `max_window_bp`.
///
/// This sorter uses `BinaryHeap<Reverse<HeapEntry>>` as a min-heap so that `peek()` and `pop()` yield
/// the smallest entry by the chosen ordering. A `BinaryHeap` is not a fully sorted container.
/// After `push`, the heap only guarantees the extremum at the top. Iterating the heap does not
/// produce sorted order. Sorted output is obtained by repeatedly popping the heap.
///
/// The algorithm buffers incoming entries and flushes anything that is safely behind a moving
/// coordinate window. Let `max_seen_start` be the largest `start` observed so far. Any entry with
/// `start <= max_seen_start - max_window_bp` cannot be preceded by future entries under the bounded
/// disorder assumption, so it can be written out immediately. At end of stream, the remaining
/// tail is at most `max_window_bp` wide and is finished with a small in-memory `Vec::sort`.
///
/// Parameters
/// ----------
/// - max_window_bp:
///     Maximum coordinate displacement the input stream can exhibit. Set this to the maximum
///     fragment length so that no future fragment starts more than `max_window_bp` bases behind
///     the current maximum start.
///
/// Behavior
/// --------
/// - Push:
///     Updates `max_seen_start`, inserts into the min-heap, then repeatedly flushes the heap
///     top while `top.start <= max_seen_start - max_window_bp`.
/// - Flush all:
///     Drains the residual heap, sorts that small tail, and writes it out to complete ordering.
/// - Ordering:
///     Writes entries strictly ordered by `(start, end, strand, qname, arrival)` within a chromosome stream.
/// - Memory bound:
///     Heap size is bounded by the number of fragments whose starts fall within the last
///     `max_window_bp` bases, not by total input size.
///
/// Complexity
/// ----------
/// - Each `push` is `O(log k)` where `k` is the current heap size within the window.
/// - Final tail sort is `O(k log k)` with `k` limited by the window, typically small.
///
/// Correctness rationale
/// ---------------------
/// - Under the bounded disorder assumption, for any observed `max_seen_start`, no unseen entry
///   can have `start < max_seen_start - max_window_bp`.
/// - Therefore any entry with `start <= max_seen_start - max_window_bp` is globally minimal among
///   the buffered entries and safe to flush without violating sorted order.
///
/// Limitations
/// -----------
/// - If the stream violates the bounded disorder assumption, output order is not guaranteed.
/// - Sorting keys are `(start, end, strand, qname, arrival)`. Extend the ordering if additional
///   keys or stability guarantees are required.
///
/// Example
/// -------
/// ```rust
/// use cfdnalab::commands::bam_to_bam::sorted_writer::{
///     RecordEntry, RecordTags, RecordWriter, WindowSorter,
/// };
/// use anyhow::Result;
/// use rust_htslib::bam::Record;
/// use std::sync::Arc;
///
/// struct Sink;
/// impl RecordWriter for Sink {
///     fn write_entry(&mut self, _entry: RecordEntry) -> Result<()> {
///         Ok(())
///     }
/// }
///
/// # fn main() -> Result<()> {
/// let mut sorter = WindowSorter::new(200);
/// let record = Record::new();
/// let entry = RecordEntry {
///     interval: cfdnalab::shared::interval::Interval::new(10, 60)?,
///     record,
///     tags: Arc::new(RecordTags {
///         fragment_length: 50,
///         coverage_weight: None,
///         fragment_count_weight: None,
///         gc_weight: None,
///     }),
/// };
/// let mut writer = Sink;
/// sorter.push(entry, &mut writer)?;
/// sorter.flush_all(&mut writer)?;
/// # Ok(())
/// # }
/// ```
///
pub struct WindowSorter {
    heap: BinaryHeap<Reverse<HeapEntry>>,
    max_window_bp: u32, // Look-back (max fragment length)
    max_seen_start: u32,
    next_serial: u64,
}

impl WindowSorter {
    /// `max_window_bp`: the maximum fragment length bound.
    pub fn new(max_window_bp: u32) -> Self {
        Self {
            heap: BinaryHeap::new(),
            max_window_bp,
            max_seen_start: 0,
            next_serial: 0,
        }
    }

    /// Push a new entry and flush anything safely behind the window to `writer`.
    pub fn push<W: RecordWriter>(&mut self, entry: RecordEntry, writer: &mut W) -> Result<()> {
        if entry.start() > self.max_seen_start {
            self.max_seen_start = entry.start();
        }
        self.heap.push(Reverse(HeapEntry {
            serial: self.next_serial,
            entry,
        }));
        self.next_serial += 1;

        let threshold = self.max_seen_start.saturating_sub(self.max_window_bp);
        while let Some(Reverse(top)) = self.heap.peek() {
            if top.start() <= threshold {
                let Reverse(e) = self.heap.pop().unwrap();
                writer.write_entry(e.entry)?;
            } else {
                break;
            }
        }
        Ok(())
    }

    /// Flush all remaining entries in sorted order.
    pub fn flush_all<W: RecordWriter>(&mut self, writer: &mut W) -> Result<()> {
        // Remaining tail is at most `max_window_bp` wide in coordinates.
        let mut rest: Vec<HeapEntry> = self.heap.drain().map(|rev| rev.0).collect();
        rest.sort(); // Sort small tail to ensure total order
        for e in rest {
            writer.write_entry(e.entry)?;
        }
        Ok(())
    }
}

/// Writer abstraction for sorted BAM output.
pub trait RecordWriter {
    fn write_entry(&mut self, entry: RecordEntry) -> Result<()>;
}

impl RecordWriter for bam::Writer {
    fn write_entry(&mut self, mut entry: RecordEntry) -> Result<()> {
        apply_aux_tags(&mut entry.record, entry.tags.as_ref())?;
        self.write(&entry.record)
            .context("writing BAM record to output")?;
        Ok(())
    }
}

fn apply_aux_tags(record: &mut Record, tags: &RecordTags) -> Result<()> {
    record
        .push_aux(FRAGMENT_LENGTH_AUX_TAG, Aux::U32(tags.fragment_length))
        .context("setting fragment_length aux tag")?;
    if let Some(weight) = tags.coverage_weight {
        record
            .push_aux(COVERAGE_WEIGHT_AUX_TAG, Aux::Float(weight))
            .context("setting coverage_weight aux tag")?;
    }
    if let Some(weight) = tags.fragment_count_weight {
        record
            .push_aux(FRAGMENT_COUNT_WEIGHT_AUX_TAG, Aux::Float(weight))
            .context("setting fragment_count_weight aux tag")?;
    }
    if let Some(weight) = tags.gc_weight {
        record
            .push_aux(GC_WEIGHT_AUX_TAG, Aux::Float(weight))
            .context("setting gc aux tag")?;
    }
    Ok(())
}

struct HeapEntry {
    serial: u64,
    entry: RecordEntry,
}

impl HeapEntry {
    #[inline]
    fn start(&self) -> u32 {
        self.entry.start()
    }

    #[inline]
    fn cmp_key(&self) -> (u32, u32, bool, &[u8], u64) {
        (
            self.entry.start(),
            self.entry.end(),
            self.entry.record.is_reverse(),
            self.entry.record.qname(),
            self.serial,
        )
    }
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_key().cmp(&other.cmp_key())
    }
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp_key() == other.cmp_key()
    }
}

impl Eq for HeapEntry {}
