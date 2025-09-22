use crate::utils::iterator_counter::FragmentCounterSnapshot;

#[derive(Debug, Default)]
pub struct GCCounters {
    /// Reads encountered
    pub total_reads: u64,
    /// Fragments collected from reads
    pub collected_fragments: u64,
    /// Forward reads accepted by first filters
    pub accepted_forward: u64,
    /// Reverse reads accepted by first filters
    pub accepted_reverse: u64,
    /// Fragments yielded from iterator
    pub yielded_fragments: u64,
    /// *Fragments* counted
    pub counted_fragments: u64,
}

impl std::ops::AddAssign for GCCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.counted_fragments += other.counted_fragments;
    }
}

impl GCCounters {
    /// Add counts from snapshot
    pub fn add_from_snapshot(&mut self, other: FragmentCounterSnapshot) {
        self.total_reads += other.incoming_reads;
        self.collected_fragments += other.produced_fragments;
        self.accepted_forward += other.accepted_forward_reads;
        self.accepted_reverse += other.accepted_reverse_reads;
        self.yielded_fragments += other.yielded_fragments;
    }
}

#[derive(Debug, Default)]
pub struct LengthsCounters {
    /// Reads encountered
    pub total_reads: u64,
    /// Fragments collected from reads
    pub collected_fragments: u64,
    /// Forward reads accepted by first filters
    pub accepted_forward: u64,
    /// Reverse reads accepted by first filters
    pub accepted_reverse: u64,
    /// Fragments yielded from iterator
    pub yielded_fragments: u64,
    /// Fragments excluded due to overlap with blacklist intervals
    pub blacklisted_fragments: u64,
    /// Fragments excluded due to extreme GC fraction
    pub gc_excl_fragments: u64,
    /// *Fragments* counted
    pub counted_fragments: u64,
}

impl std::ops::AddAssign for LengthsCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.yielded_fragments += other.yielded_fragments;
        self.blacklisted_fragments += other.blacklisted_fragments;
        self.gc_excl_fragments += other.gc_excl_fragments;
        self.counted_fragments += other.counted_fragments;
    }
}

impl LengthsCounters {
    /// Add counts from snapshot
    pub fn add_from_snapshot(&mut self, other: FragmentCounterSnapshot) {
        self.total_reads += other.incoming_reads;
        self.collected_fragments += other.produced_fragments;
        self.accepted_forward += other.accepted_forward_reads;
        self.accepted_reverse += other.accepted_reverse_reads;
        self.yielded_fragments += other.yielded_fragments;
    }
}

#[derive(Debug, Default)]
pub struct NormalizeGenomeCounters {
    /// Reads encountered
    pub total_reads: u64,
    /// Fragments collected from reads
    pub collected_fragments: u64,
    /// Forward reads accepted by first filters
    pub accepted_forward: u64,
    /// Reverse reads accepted by first filters
    pub accepted_reverse: u64,
    /// Fragments yielded from iterator
    pub yielded_fragments: u64,
    /// *Fragments* counted
    pub counted_fragments: u64,
}

impl std::ops::AddAssign for NormalizeGenomeCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.yielded_fragments += other.yielded_fragments;
        self.counted_fragments += other.counted_fragments;
    }
}

impl NormalizeGenomeCounters {
    /// Add counts from snapshot
    pub fn add_from_snapshot(&mut self, other: FragmentCounterSnapshot) {
        self.total_reads += other.incoming_reads;
        self.collected_fragments += other.produced_fragments;
        self.accepted_forward += other.accepted_forward_reads;
        self.accepted_reverse += other.accepted_reverse_reads;
        self.yielded_fragments += other.yielded_fragments;
    }
}

#[derive(Debug, Default)]
pub struct FCoverageCounters {
    /// Reads encountered
    pub total_reads: u64,
    /// Fragments collected from reads
    pub collected_fragments: u64,
    /// Forward reads accepted by first filters
    pub accepted_forward: u64,
    /// Reverse reads accepted by first filters
    pub accepted_reverse: u64,
    /// *Fragments* counted
    pub counted_fragments: u64,
}

impl std::ops::AddAssign for FCoverageCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.counted_fragments += other.counted_fragments;
    }
}

impl FCoverageCounters {
    /// Add counts from snapshot
    pub fn add_from_snapshot(&mut self, other: FragmentCounterSnapshot) {
        self.total_reads += other.incoming_reads;
        self.collected_fragments += other.produced_fragments;
        self.accepted_forward += other.accepted_forward_reads;
        self.accepted_reverse += other.accepted_reverse_reads;
        self.counted_fragments += other.yielded_fragments;
    }
}

/* Profile Groups */

#[derive(Debug, Default)]
pub struct ProfileGroupsCounters {
    /// Reads encountered
    pub total_reads: u64,
    /// Fragments collected from reads
    pub collected_fragments: u64,
    /// Forward reads accepted by first filters
    pub accepted_forward: u64,
    /// Reverse reads accepted by first filters
    pub accepted_reverse: u64,
    /// Fragments excluded due to overlap with blacklist intervals
    pub blacklisted_fragments: u64,
    /// *Fragments* created (double counts likely)
    pub yielded_fragments: u64,
    /// *Fragments* counted
    pub counted_fragments: u64,
}

impl std::ops::AddAssign for ProfileGroupsCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.blacklisted_fragments += other.blacklisted_fragments;
        self.yielded_fragments += other.yielded_fragments;
        self.counted_fragments += other.counted_fragments;
    }
}

impl ProfileGroupsCounters {
    /// Add counts from snapshot
    pub fn add_from_snapshot(&mut self, other: FragmentCounterSnapshot) {
        self.total_reads += other.incoming_reads;
        self.collected_fragments += other.produced_fragments;
        self.accepted_forward += other.accepted_forward_reads;
        self.accepted_reverse += other.accepted_reverse_reads;
        self.yielded_fragments += other.yielded_fragments;
    }
}
