#[cfg(any(feature = "cmd_ends", feature = "cmd_lengths"))]
pub(crate) mod cigar_counts;
#[cfg(feature = "cmd_ends")]
pub(crate) mod ends_fragment;
#[cfg(feature = "cmd_bam_to_frag")]
pub(crate) mod frag_file_fragment;
#[cfg(feature = "cmd_lengths")]
pub(crate) mod indel_counting_fragment;
pub(crate) mod minimal_fragment;
#[cfg(feature = "cmd_fcoverage")]
pub(crate) mod segment_fragment;
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) mod segment_kmer_fragment;
#[cfg(feature = "cmd_bam_to_bam")]
pub(crate) mod with_records_fragment;
