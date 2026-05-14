use crate::{
    commands::cli_common::{
        ApplyGCArgs, ChromosomeArgs, IOCArgs, LoggingArgs, ScaleGenomeArgs, UnpairedArgs,
        resolve_length_bin_edges,
    },
    commands::midpoints::smoothing::MidpointSmoothing,
    shared::{
        blacklist::BlacklistStrategy,
        constants::{MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION},
    },
};
use anyhow::Result;
use std::path::PathBuf;

/// Count positional fragment **midpoint** coverage in groups of genomic windows.
///
/// **Midpoints**: The center of the fragment span, with ties (in even-sized fragments)
/// randomly and reproducibly assigned to the left or right mid-position to avoid bias
/// from always rounding in the same direction.
///
/// **Groups**: The coverage profiles are collapsed (summed per position) across windows in a group.
/// E.g., transcription factors as groups with binding sites as windows, yielding the
/// overall midpoint profile per transcription factor.
///
/// **Smoothing**: Final profiles are smoothed with an order-3 Savitzky-Golay
/// filter unless `--smoothing none` is used.
///
/// **Strandedness**: When the `--intervals` carry strand information (`+`/`-`/`.`),
/// reverse-stranded (`-`) intervals write into group profiles in reverse order.
///
/// ## Fragment span definition
///
/// **Paired-end**: `[forward.pos, reverse.reference_end)`, the reference span
/// from the first aligned position on the forward read to the last aligned
/// position on the reverse read.
///
/// **Unpaired** where each read is a fragment: `[read.pos, read.reference_end)`,
/// the reference span from the first to the last aligned position on the read.
///
/// The utilized fragment length range is specified via `--length-bins`.
///
/// ## Always-on exclusion criteria
///
/// The following criteria always exclude a read:
///
/// The read is secondary, supplementary or duplicate.
/// The read failed quality check.
///
/// **Paired-end input only**:
/// The read or mate read is unmapped.
/// The read is mapped to a different `tid` than the mate.
/// The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct MidpointsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls of the command.
    ///
    /// Examples produce files like:
    ///   `<prefix>.midpoint_profiles.npy`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// The grouped fixed-size intervals to count within `[path]`
    ///
    /// A BED-like file of genomic intervals, their respective group names, and optionally
    /// the interval strandedness.
    ///
    /// Must be sorted by the `chromosome` and `start` coordinates, and
    /// all intervals must have the same size.
    ///
    /// Sites with the same group name are collapsed to a single profile.
    ///
    /// Strand tokens are `+`, `-`, or `.`. With six or more columns, only column 6 is read as
    /// strand. With exactly five columns, column 5 may be strand. If column 5 looks stranded but
    /// column 6 exists and does not, the file is rejected as ambiguous.
    ///
    /// Forward and unknown strandedness use normal genomic order while reverse stranded intervals
    /// write into the grouped profiles in reverse order.
    ///
    /// Required columns: `chromosome, start, end, group_name`. No header.
    /// Optional fifth and sixth columns follow the strand rules above.
    ///
    /// Note: Besides chromosome-filtering and blacklist-filtering (see explanation in
    /// `--keep-blacklisted-intervals`), no additional interval filtering is performed.
    /// It is up to the user to remove duplicate intervals, within-group overlapping intervals,
    /// and other potential clean-ups beforehand.
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'w',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub intervals: PathBuf,

    /// Edges of fragment length bins to count in `[string(s)]`
    ///
    /// This also defines the min and max fragment lengths.
    ///
    /// Accepted forms:
    ///
    /// - A single value with `start:end:step`:
    ///   Creates contiguous bins from `start` to `end` (end-exclusive) in `step` increments.
    ///   Example: `30:1000:10` -> bins `[30,40), [40,50), ..., [990,1000)`.
    ///
    /// - Multiple integer values interpreted as bin edges:
    ///   Example: `--length-bins 30 80 150 220 500 1001` -> bins `[30,80), [80,150), ..., [500,1001)`.
    ///
    /// **NOTE**: Memory consumption increases linearly with the number of bins.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            default_values_t = [String::from("30"), String::from("1001")],
            help_heading = "Core"
        )
    )]
    pub length_bins: Vec<String>,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    /// Average adjacent positions in bins after counting and smoothing `[integer]`
    ///
    /// Defaults to 1 for full resolution.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1", value_parser = clap::value_parser!(u32).range(1..), help_heading = "Post-processing")
    )]
    pub bin_size: u32,

    /// Smooth final midpoint profiles with an order-3 Savitzky-Golay filter `[string]`
    ///
    /// By default smoothing is applied with a 165 bp window (nucleosome-scale).
    /// This was selected based on the defaults in Griffin.
    ///
    /// Use `none` to write unsmoothed profiles.
    ///
    /// Use `savgol=<odd_bp>` to set a different odd window size in base pairs. For example,
    /// `--smoothing savgol=155` uses a 155 bp smoothing window.
    ///
    /// When smoothing is active, intervals are counted with additional flank positions to avoid
    /// edge effects in the smoothed values. Flanked intervals cannot exceed chromosome boundaries.
    /// After smoothing final grouped profiles, the command trims the flanks and writes exactly
    /// the positions from `--intervals`.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "savgol=165", help_heading = "Post-processing")
    )]
    pub smoothing: MidpointSmoothing,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    ///
    /// This is **NOT** recommended by default as it trims the tails of the length distribution.
    ///
    /// Note, that we only keep inward-directed fragments within the specified length range, so
    /// there's no real need for proper-pair filtering.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// **NOTE**: It may be an advantage to instead remove intervals that lie within
    /// half the maximum fragment length of blacklisted regions from the `--intervals` file.
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-min-size",
            default_value = "1",
            help_heading = "Filtering"
        )
    )]
    pub blacklist_min_size: u64,

    /// The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
    ///
    /// Possible values:
    ///     `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
    ///
    /// `midpoint` checks the single central base for odd fragments and either
    /// central base for even fragments.
    ///
    /// Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-strategy",
            default_value = "any",
            ignore_case = true,
            help_heading = "Filtering"
        )
    )]
    pub blacklist_strategy: BlacklistStrategy,

    /// Don't filter out intervals that overlap nearby blacklisted regions `[flag]`
    ///
    /// **Edge bias**: Fragments overlapping blacklisted bases are filtered before
    /// they can contribute midpoint counts. This can create artificial
    /// dips near profile edges if an interval is close enough to a blacklist for relevant
    /// fragments to be removed on one side but not the other.
    ///
    /// To avoid that edge bias, intervals within `ceil(max_fragment_length / 2) + smoothing_flank`
    /// from blacklisted regions are removed prior to counting.
    /// `smoothing_flank` is half the Savitzky-Golay window when smoothing is active and
    /// `0` when `--smoothing none`.
    ///
    /// Set this flag to keep those intervals in the site set.
    /// Fragment-level blacklist filtering still applies.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub keep_blacklisted_intervals: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub gc: ApplyGCArgs,

    /// Optional 2bit reference genome file [path]
    ///
    /// NOTE: Required for GC correction, otherwise ignored.
    ///
    /// E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'r',
            long,
            value_parser,
            required = false,
            help_heading = "GC Correction"
        )
    )]
    pub ref_2bit: Option<PathBuf>,

    /// Group indices to plot as midpoint profiles `[integers]`
    ///
    /// Comma separated list of zero-based group indices to plot after counting.
    ///
    /// This plotting step is intended for quick QC of the outputs. It's not
    /// optimized for publication etc. (although feel free!)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_delimiter = ',',
            default_values_t = [0_usize],
            help_heading = "Plotting"
        )
    )]
    pub plot_groups: Vec<usize>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl MidpointsConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs, intervals: PathBuf) -> Self {
        Self {
            ioc,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix: String::new(),
            intervals,
            length_bins: vec!["30".to_string(), "1001".to_string()],
            bin_size: 1,
            smoothing: MidpointSmoothing::default(),
            tile_size: 20000000,
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
            keep_blacklisted_intervals: false,
            gc: ApplyGCArgs {
                gc_file: None,
                gc_tag: None,
                neutralize_invalid_gc: false,
            },
            ref_2bit: None,
            plot_groups: vec![0],
            logging: LoggingArgs::default(),
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_length_bins(&mut self, edges: Vec<u32>) {
        assert!(
            edges.len() >= 2,
            "length bin edges must contain at least two values"
        );
        self.length_bins = edges.into_iter().map(|edge| edge.to_string()).collect();
    }

    pub fn set_length_bins_spec<S: Into<String>>(&mut self, spec: S) {
        self.length_bins = vec![spec.into()];
    }

    pub fn set_bin_size(&mut self, bin_size: u32) {
        assert!(bin_size >= 1, "bin_size must be at least 1");
        self.bin_size = bin_size;
    }

    pub fn set_smoothing(&mut self, smoothing: MidpointSmoothing) {
        self.smoothing = smoothing;
    }

    pub fn resolve_length_bins(&self) -> Result<Vec<u32>> {
        resolve_length_bin_edges(
            &self.length_bins,
            MIN_ACGT_BASES_FOR_GC_FRACTION,
            MAX_SUPPORTED_FRAGMENT_LENGTH,
        )
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_keep_blacklisted_intervals(&mut self, keep: bool) {
        self.keep_blacklisted_intervals = keep;
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.scale_genome = scale;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}

#[cfg(test)]
mod tests {
    include!("config_tests.rs");
}
