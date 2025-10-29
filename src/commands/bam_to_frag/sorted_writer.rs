use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::io::{self, Write};

/// A buffered fragment line with its sort keys.
#[derive(Eq, PartialEq)]
pub struct Entry {
    pub start: u32,
    pub end: u32,
    pub line: String,
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Ascending by (start, end, line). We use Reverse<Entry> to make a min-heap.
        (self.start, self.end, &self.line).cmp(&(other.start, other.end, &other.line))
    }
}
impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Maintains a bounded in-memory buffer to emit fragments in strict (start, end, line) order
/// from an almost-sorted stream whose disorder is limited by `max_window_bp`.
///
/// This sorter uses `BinaryHeap<Reverse<Entry>>` as a min-heap so that `peek()` and `pop()` yield
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
///     Emits entries strictly ordered by `(start, end, line)` within a chromosome stream.
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
/// - Sorting keys are `(start, end, line)`. Extend the `Entry` ordering if additional keys
///   or stability guarantees are required.
///
/// Example
/// -------
/// ```rust
/// let mut sorter = WindowSorter::new(max_fragment_length);
/// // Stream fragments:
/// sorter.push(Entry { start, end, line }, &mut writer)?;
/// // End of stream:
/// sorter.flush_all(&mut writer)?;
/// ```

pub struct WindowSorter {
    pub heap: BinaryHeap<Reverse<Entry>>,
    pub max_window_bp: u32, // Look-back (max fragment length)
    pub max_seen_start: u32,
}

impl WindowSorter {
    /// `max_window_bp`: the maximum fragment length bound.
    pub fn new(max_window_bp: u32) -> Self {
        Self {
            heap: BinaryHeap::new(),
            max_window_bp,
            max_seen_start: 0,
        }
    }

    /// Push a new entry and flush anything safely behind the window to `w`.
    pub fn push<W: Write>(&mut self, entry: Entry, w: &mut W) -> io::Result<()> {
        if entry.start > self.max_seen_start {
            self.max_seen_start = entry.start;
        }
        self.heap.push(Reverse(entry));

        let threshold = self.max_seen_start.saturating_sub(self.max_window_bp);
        while let Some(Reverse(top)) = self.heap.peek() {
            if top.start <= threshold {
                let Reverse(e) = self.heap.pop().unwrap();
                w.write_all(e.line.as_bytes())?;
            } else {
                break;
            }
        }
        Ok(())
    }

    /// Flush all remaining entries in sorted order.
    pub fn flush_all<W: Write>(&mut self, w: &mut W) -> io::Result<()> {
        // Remaining tail is at most `max_window_bp` wide in coordinates.
        let mut rest: Vec<Entry> = self.heap.drain().map(|rev| rev.0).collect();
        rest.sort(); // Sort small tail to ensure total order
        for e in rest {
            w.write_all(e.line.as_bytes())?;
        }
        Ok(())
    }
}
