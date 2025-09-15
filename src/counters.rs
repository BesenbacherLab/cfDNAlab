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
    /// Fragments excluded due to overlap with blacklist intervals
    pub blacklisted_fragments: u64,
    /// Fragments excluded due to extreme GC fraction
    pub gc_excl_fragments: u64,
    // Fragments excluded for being too short/long
    pub illegal_length_fragments: u64,
    /// *Fragments* counted
    pub counted_fragments: u64,
}

impl std::ops::AddAssign for LengthsCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.blacklisted_fragments += other.blacklisted_fragments;
        self.illegal_length_fragments += other.illegal_length_fragments;
        self.gc_excl_fragments += other.gc_excl_fragments;
        self.counted_fragments += other.counted_fragments;
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
    // Fragments excluded for being too short/long
    pub illegal_length_fragments: u64,
    /// *Fragments* counted
    pub counted_fragments: u64,
}

impl std::ops::AddAssign for NormalizeGenomeCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.illegal_length_fragments += other.illegal_length_fragments;
        self.counted_fragments += other.counted_fragments;
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

// pub struct FragmentCounterSnapshot {
//     pub incoming_reads: u64,
//     pub incoming_fragments: u64,
//     pub accepted_forward_reads: u64,
//     pub accepted_reverse_reads: u64,
//     pub produced_fragments: u64,
//     pub yielded_fragments: u64,
// }
