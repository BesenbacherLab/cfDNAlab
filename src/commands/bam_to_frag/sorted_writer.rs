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

/// Sorts an almost-sorted fragment stream using a bounded coordinate window.
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
