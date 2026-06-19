#![allow(
    dead_code,
    reason = "feature-limited builds compile shared counter definitions even when no enabled command constructs them"
)]

use crate::shared::iterator_counter::FragmentCounterSnapshot;
use std::ops::AddAssign;

#[derive(Debug, Default, Clone, Copy)]
pub struct BaseCounters {
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

impl AddAssign for BaseCounters {
    fn add_assign(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.collected_fragments += other.collected_fragments;
        self.accepted_forward += other.accepted_forward;
        self.accepted_reverse += other.accepted_reverse;
        self.yielded_fragments += other.yielded_fragments;
        self.counted_fragments += other.counted_fragments;
    }
}

impl BaseCounters {
    /// Add counts from snapshot
    pub(crate) fn add_from_snapshot(&mut self, snap: FragmentCounterSnapshot) {
        self.total_reads += snap.incoming_reads;
        self.collected_fragments += snap.produced_fragments;
        self.accepted_forward += snap.accepted_forward_reads;
        self.accepted_reverse += snap.accepted_reverse_reads;
        self.yielded_fragments += snap.yielded_fragments;
    }
}

/// Macro to declare a counters struct with a BaseCounters + extra fields,
/// plus AddAssign and add_from_snapshot impls.
#[allow(
    unused_macros,
    reason = "feature-limited builds may compile this helper without any enabled counter declarations"
)]
macro_rules! counter_struct {
    // No extra fields
    ($name:ident ;) => {
        #[derive(Debug, Default, Clone, Copy)]
        pub struct $name {
            pub base: BaseCounters,
        }
        impl AddAssign for $name {
            fn add_assign(&mut self, other: Self) {
                self.base += other.base;
            }
        }
        impl $name {
            /// Add counts from snapshot
            #[allow(dead_code)]
            pub(crate) fn add_from_snapshot(&mut self, snap: FragmentCounterSnapshot) {
                self.base.add_from_snapshot(snap);
            }
        }
    };
    // With extra fields
    ($name:ident ; $( $field:ident : $ty:ty ),+ $(,)? ) => {
        #[derive(Debug, Default, Clone, Copy)]
        pub struct $name {
            pub base: BaseCounters,
            $( pub $field: $ty, )+
        }
        impl AddAssign for $name {
            fn add_assign(&mut self, other: Self) {
                self.base += other.base;
                $( self.$field += other.$field; )+
            }
        }
        impl $name {
            /// Add counts from snapshot
            #[allow(dead_code)]
            pub(crate) fn add_from_snapshot(&mut self, snap: FragmentCounterSnapshot) {
                self.base.add_from_snapshot(snap);
            }
        }
    };
}

// Declarations

#[cfg(feature = "cmd_gc_bias")]
counter_struct!(GCCounters;);

#[cfg(feature = "cmd_fragment_kmers")]
counter_struct!(FragmentKmersCounters;
    blacklisted_fragments: u64,
    gc_failed_fragments: u64,
    gc_out_of_range_tags: u64,
    gc_missing_tags: u64
);

#[cfg(feature = "cmd_fcoverage")]
counter_struct!(FCoverageCounters;
    tile_owned_normalization_fragments: u64,
    tile_owned_normalization_length_sum: u64,
    gc_failed_fragments: u64,
    gc_out_of_range_tags: u64,
    gc_missing_tags: u64
);

#[cfg(feature = "cmd_wps")]
counter_struct!(WPSCounters;
    gc_failed_fragments: u64,
    gc_out_of_range_tags: u64,
    gc_missing_tags: u64
);

#[cfg(feature = "cmd_wps_peaks")]
counter_struct!(WPSPeaksCounters;
    gc_failed_fragments: u64,
    gc_out_of_range_tags: u64,
    gc_missing_tags: u64
);

#[cfg(feature = "cmd_wps_peaks")]
impl From<WPSCounters> for WPSPeaksCounters {
    fn from(other: WPSCounters) -> Self {
        Self {
            base: other.base,
            gc_failed_fragments: other.gc_failed_fragments,
            gc_out_of_range_tags: other.gc_out_of_range_tags,
            gc_missing_tags: other.gc_missing_tags,
        }
    }
}

#[cfg(feature = "cmd_lengths")]
counter_struct!(LengthsCounters;
    blacklisted_fragments: u64,
    gc_failed_fragments: u64
);

#[cfg(feature = "cmd_ends")]
counter_struct!(EndsCounters;
    blacklisted_fragments: u64,
    gc_failed_fragments: u64,
    gc_out_of_range_tags: u64,
    gc_missing_tags: u64,
    counted_motifs: u64
);

#[cfg(feature = "cmd_midpoints")]
counter_struct!(ProfileGroupsCounters;
    blacklisted_fragments: u64,
    gc_failed_fragments: u64,
    gc_out_of_range_tags: u64,
    gc_missing_tags: u64,
);

#[cfg(feature = "cmd_bam_to_bam")]
counter_struct!(BamToBamCounters; blacklisted_fragments: u64, gc_failed_fragments: u64);

#[cfg(feature = "cmd_bam_to_frag")]
counter_struct!(BamToFragCounters; blacklisted_fragments: u64, gc_failed_fragments: u64);
