pub(crate) mod bam;
pub(crate) mod base;
pub(crate) mod bed;
pub(crate) mod blacklist;
pub(crate) mod cli_output;
pub(crate) mod clip_mode;
pub(crate) mod constants;
#[cfg(feature = "cmd_fcoverage")]
pub(crate) mod coverage;
#[cfg(any(feature = "cmd_lengths", feature = "cmd_fcoverage"))]
pub(crate) mod formatters;
pub(crate) mod fragment;
#[cfg(any(
    feature = "cmd_bam_to_bam",
    feature = "cmd_bam_to_frag",
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_gc_bias",
    feature = "cmd_lengths",
    feature = "cmd_midpoints",
    feature = "cmd_wps",
    feature = "cmd_wps_peaks"
))]
pub(crate) mod fragment_iterators;
pub(crate) mod gc_tag;
pub(crate) mod indel_mode;
pub(crate) mod interval;
pub(crate) mod io;
pub(crate) mod iterator_counter;
#[cfg(any(feature = "cmd_ends", feature = "cmd_fragment_kmers"))]
pub(crate) mod kmers;
#[cfg(any(feature = "cmd_lengths", feature = "cmd_midpoints"))]
pub(crate) mod length_axis;
pub(crate) mod logging;
#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_gc_bias",
    feature = "cmd_lengths",
    feature = "cmd_midpoints"
))]
pub(crate) mod midpoint;
pub(crate) mod overlaps;
pub(crate) mod positioning;
#[cfg(uses_progress_reporting)]
pub(crate) mod progress;
pub(crate) mod read;
pub(crate) mod reference;
#[cfg(feature = "cmd_gc_bias")]
pub(crate) mod sampling;
pub(crate) mod scale_genome;
#[cfg(uses_temp_chrom_names)]
pub(crate) mod temp_chrom_names;
pub(crate) mod thread_pool;
pub(crate) mod tiled_run;
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) mod visualization;
#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths",
    feature = "cmd_midpoints"
))]
pub(crate) mod window_fetch;
#[cfg(any(
    feature = "cmd_bam_to_bam",
    feature = "cmd_bam_to_frag",
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths",
    feature = "cmd_gc_bias",
    feature = "cmd_wps",
    feature = "cmd_wps_peaks"
))]
pub(crate) mod windowing;
#[cfg(any(feature = "cmd_bam_to_frag", feature = "cmd_fcoverage"))]
pub(crate) mod writers;
#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_gc_bias",
    feature = "cmd_midpoints"
))]
pub(crate) mod zarr;
// Plotting helpers gated behind the plotters feature
#[cfg(feature = "plotters")]
pub(crate) mod plotters;
