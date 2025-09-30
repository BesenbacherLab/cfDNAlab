// Read + fragment counters for the generic fragment iterator
// Local: For single-thread use
// Shared: Safe for multi-thread use

#[derive(Debug, Clone, Copy, Default)]
pub struct FragmentCounterSnapshot {
    pub incoming_reads: u64,
    pub incoming_fragments: u64,
    pub accepted_forward_reads: u64,
    pub accepted_reverse_reads: u64,
    pub produced_fragments: u64,
    pub yielded_fragments: u64,
}

// Generic counter trait
pub trait FragmentCounters: Send + 'static {
    fn inc_incoming_reads(&mut self);
    fn inc_incoming_fragments(&mut self);
    fn inc_accepted_reads(&mut self, is_reverse: bool);
    fn inc_produced_fragments(&mut self);
    fn inc_yielded_fragments(&mut self);
    fn snapshot(&self) -> FragmentCounterSnapshot;
}

/// Zero-overhead default for no-stats runs
#[derive(Clone, Default)]
pub struct NoopCounters;

impl FragmentCounters for NoopCounters {
    #[inline]
    fn inc_incoming_reads(&mut self) {}
    #[inline]
    fn inc_incoming_fragments(&mut self) {}
    #[inline]
    fn inc_accepted_reads(&mut self, _is_reverse: bool) {}
    #[inline]
    fn inc_produced_fragments(&mut self) {}
    #[inline]
    fn inc_yielded_fragments(&mut self) {}
    #[inline]
    fn snapshot(&self) -> FragmentCounterSnapshot {
        FragmentCounterSnapshot::default()
    }
}

/* Single-thread */

/// Fast, single-threaded counters for single-thread iterators
#[derive(Clone, Default)]
pub struct LocalCounters {
    incoming_reads: u64,
    incoming_fragments: u64,
    accepted_forward_reads: u64,
    accepted_reverse_reads: u64,
    produced_fragments: u64,
    yielded_fragments: u64,
}

impl LocalCounters {
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }
}

impl FragmentCounters for LocalCounters {
    #[inline]
    fn inc_incoming_reads(&mut self) {
        self.incoming_reads += 1;
    }
    #[inline]
    fn inc_incoming_fragments(&mut self) {
        self.incoming_fragments += 1;
    }
    #[inline]
    fn inc_accepted_reads(&mut self, is_reverse: bool) {
        if is_reverse {
            self.accepted_reverse_reads += 1;
        } else {
            self.accepted_forward_reads += 1;
        }
    }
    #[inline]
    fn inc_produced_fragments(&mut self) {
        self.produced_fragments += 1;
    }
    #[inline]
    fn inc_yielded_fragments(&mut self) {
        self.yielded_fragments += 1;
    }
    #[inline]
    fn snapshot(&self) -> FragmentCounterSnapshot {
        FragmentCounterSnapshot {
            incoming_reads: self.incoming_reads,
            incoming_fragments: self.incoming_fragments,
            accepted_forward_reads: self.accepted_forward_reads,
            accepted_reverse_reads: self.accepted_reverse_reads,
            produced_fragments: self.produced_fragments,
            yielded_fragments: self.yielded_fragments,
        }
    }
}

/* Cross-thread */

/// Cross-thread, shareable counters (atomics). Use when multiple iterators should report to one place.
#[derive(Clone, Default)]
pub struct SharedCounters {
    inner: std::sync::Arc<SharedCountersInner>,
}

#[derive(Default)]
struct SharedCountersInner {
    incoming_reads: std::sync::atomic::AtomicU64,
    incoming_fragments: std::sync::atomic::AtomicU64,
    accepted_forward_reads: std::sync::atomic::AtomicU64,
    accepted_reverse_reads: std::sync::atomic::AtomicU64,
    produced_fragments: std::sync::atomic::AtomicU64,
    yielded_fragments: std::sync::atomic::AtomicU64,
}

impl SharedCounters {
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }
}

impl FragmentCounters for SharedCounters {
    #[inline]
    fn inc_incoming_reads(&mut self) {
        use std::sync::atomic::Ordering::Relaxed;
        self.inner.incoming_reads.fetch_add(1, Relaxed);
    }
    #[inline]
    fn inc_incoming_fragments(&mut self) {
        use std::sync::atomic::Ordering::Relaxed;
        self.inner.incoming_fragments.fetch_add(1, Relaxed);
    }
    #[inline]
    fn inc_accepted_reads(&mut self, is_reverse: bool) {
        use std::sync::atomic::Ordering::Relaxed;
        if is_reverse {
            self.inner.accepted_reverse_reads.fetch_add(1, Relaxed);
        } else {
            self.inner.accepted_forward_reads.fetch_add(1, Relaxed);
        }
    }
    #[inline]
    fn inc_produced_fragments(&mut self) {
        use std::sync::atomic::Ordering::Relaxed;
        self.inner.produced_fragments.fetch_add(1, Relaxed);
    }
    #[inline]
    fn inc_yielded_fragments(&mut self) {
        use std::sync::atomic::Ordering::Relaxed;
        self.inner.yielded_fragments.fetch_add(1, Relaxed);
    }
    #[inline]
    fn snapshot(&self) -> FragmentCounterSnapshot {
        use std::sync::atomic::Ordering::Relaxed;
        FragmentCounterSnapshot {
            incoming_reads: self.inner.incoming_reads.load(Relaxed),
            incoming_fragments: self.inner.incoming_fragments.load(Relaxed),
            accepted_forward_reads: self.inner.accepted_forward_reads.load(Relaxed),
            accepted_reverse_reads: self.inner.accepted_reverse_reads.load(Relaxed),
            produced_fragments: self.inner.produced_fragments.load(Relaxed),
            yielded_fragments: self.inner.yielded_fragments.load(Relaxed),
        }
    }
}
