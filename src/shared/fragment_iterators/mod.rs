// This module unifies BAM-based pairing and ready-made fragments under
// a single iterator interface, split by fragment payload type so the wiring stays readable.
//
// Usage examples live in `core.rs`, where the shared `PairingAdapter` and counter plumbing live.

mod basic;
mod core;
mod ends;
mod frag_file;
mod indel_counts;
mod kmer_segments;
mod segments;
mod with_records;

pub use basic::*;
pub use core::*;
pub use ends::*;
pub use frag_file::*;
pub use indel_counts::*;
pub use kmer_segments::*;
pub use segments::*;
pub use with_records::*;
