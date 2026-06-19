//! Test helpers for building small cfDNAlab input files and reading public outputs.
//!
//! This module is for downstream crates that want to test code which calls
//! cfDNAlab APIs or cfDNAlab command-style runners. Enable it with the
//! `testing` cargo feature:
//!
//! ```toml
//! cfdnalab = { version = "...", features = ["testing"] }
//! ```
//!
//! The builders create real files in temporary directories, for example BAM,
//! BAI, two-bit, BED, and TSV inputs. Returned `Temp*` values own those
//! directories, so paths remain valid while the value is alive and are removed
//! when it is dropped.
//!
//! The helpers are meant to make expected values derivable at the test site.
//! They use small explicit coordinates, validate fragment spans, and document
//! assumptions such as contig names, fragment spans, CIGAR operations, and
//! sequence content. They are not a production API for running analyses.

pub mod bam;
pub mod bed;
pub mod gc_packages;
pub mod output_readers;
pub mod reference;
pub mod scaling;

pub use bam::{
    Cigar, FragmentSpec, PairedFragmentSpec, ReadNamePolicy, ReadSpec, TempBam, TempBamBuilder,
    bam_from_fragment_starts, bam_from_fragments, bam_from_fragments_with_record_indexed_names,
    bam_with_indel_and_softclip_reads, long_inward_fragment_series_bam,
    single_contig_inward_pair_bam, single_read_bam_with_qualities,
};
pub use bed::{Bed4Row, write_bed4};
#[cfg(feature = "cmd_gc_bias")]
pub use gc_packages::{
    build_command_produced_gc_correction_package_for_length,
    build_command_produced_gc_correction_package_for_range,
    build_command_produced_gc_correction_package_from_reference_windows,
    build_command_produced_gc_correction_package_from_reference_windows_for_range,
};
#[cfg(feature = "cmd_gc_bias")]
pub use gc_packages::{
    write_constant_gc_correction_package, write_two_bin_gc_correction_package,
    write_unit_gc_correction_package, write_unit_gc_correction_package_for_range,
};
pub use output_readers::{
    ReferenceGCPackageMetadata, ReferenceGCPackageOutput, read_length_counts_text,
    read_length_counts_tsv, read_midpoint_zarr_counts, read_midpoint_zarr_i32_1d,
    read_midpoint_zarr_u32_1d, read_reference_gc_package, read_zst_to_string, touch_file,
};
pub use reference::{
    RepeatingContigSpec, TempTwoBit, twobit_from_sequences, twobit_with_single_repeating_contig,
};
pub use scaling::{ScalingFactorRow, write_scaling_factors_tsv};
