#[cfg(any(feature = "cmd_ends", feature = "cmd_fragment_kmers"))]
pub(crate) mod kmer_codec;
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) mod nearest_guard;
#[cfg(any(feature = "cmd_ends", feature = "cmd_fragment_kmers"))]
pub(crate) mod process_counts;
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) mod write;
