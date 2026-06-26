#[cfg(uses_kmers)]
pub(crate) mod kmer_codec;
#[cfg(any(feature = "cmd_ends", feature = "cmd_ref_kmers"))]
pub(crate) mod motifs_file;
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) mod nearest_guard;
#[cfg(uses_kmers)]
pub(crate) mod process_counts;
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) mod write;
