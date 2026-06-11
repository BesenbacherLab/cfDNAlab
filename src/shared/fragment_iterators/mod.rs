// This module unifies BAM-based pairing and ready-made fragments under
// a single iterator interface, split by fragment record type so the wiring stays readable.
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
#[cfg(feature = "cmd_fragment_kmers")]
mod kmer_segments;
#[cfg(feature = "cmd_allelic_fragments")]
mod read_order;
mod segments;
#[cfg(any(feature = "cmd_bam_to_bam", feature = "cmd_allelic_fragments"))]
mod with_records;

pub(crate) use basic::*;
pub(crate) use core::*;
#[cfg(feature = "cmd_ends")]
pub(crate) use ends::*;
#[cfg(feature = "cmd_bam_to_frag")]
pub(crate) use frag_file::*;
#[cfg(feature = "cmd_lengths")]
pub(crate) use indel_counts::*;
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) use kmer_segments::*;
#[cfg(feature = "cmd_allelic_fragments")]
pub(crate) use read_order::*;
pub(crate) use segments::*;
#[cfg(any(feature = "cmd_bam_to_bam", feature = "cmd_allelic_fragments"))]
pub(crate) use with_records::*;
