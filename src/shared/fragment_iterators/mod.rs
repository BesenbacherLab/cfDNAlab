// This module unifies BAM-based pairing and ready-made fragments under
// a single iterator interface, split by fragment payload type so the wiring stays readable.
//
// Usage examples live in `core.rs`, where the shared `PairingAdapter` and counter plumbing live.

mod basic;
mod core;
#[cfg(feature = "cmd_ends")]
mod ends;
#[cfg(feature = "cmd_bam_to_frag")]
mod frag_file;
#[cfg(feature = "cmd_lengths")]
mod indel_counts;
mod kmer_segments;
mod segments;
#[cfg(feature = "cmd_bam_to_bam")]
mod with_records;

pub use basic::*;
pub use core::*;
#[cfg(feature = "cmd_ends")]
pub use ends::*;
#[cfg(feature = "cmd_bam_to_frag")]
pub use frag_file::*;
#[cfg(feature = "cmd_lengths")]
pub use indel_counts::*;
pub use kmer_segments::*;
pub use segments::*;
#[cfg(feature = "cmd_bam_to_bam")]
pub use with_records::*;
