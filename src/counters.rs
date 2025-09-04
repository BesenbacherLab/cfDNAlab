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
