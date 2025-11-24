use std::path::PathBuf;

use crate::commands::cli_common::{ApplyGCArgs, ScaleGenomeArgs};
use crate::commands::cli_common::{ChromosomeArgs, FragmentLengthArgs, IOCArgs, WindowsArgs};
use crate::commands::fcoverage::window_results::CoverageWindowAction;

/// Count positional **fragment** coverage across the genome.
///
/// Only paired-end fragments with both reads present are counted. By default,
/// the entire fragment span `[start(forward), end(reverse))` is counted, except for
/// deletions and skipped regions that are not covered by the other read.
///
/// ## Windowing
///
/// When specifying windows (`--by-bed` or `--by-size`), one of the following outputs
/// is possible:
///
///  - Get the average (default) or total coverage per window.
///
///  - Get the positional coverage for the included windows only (`--by-bed` *only*).
///    Excludes all positions that do not overlap a window from the output.
///    Choose between:
///     1) Indexed: Adds the original window index as an output column and keeps duplicate positions.
///     2) Unique: Overlapping windows are merged to avoid duplicate positions.
///
/// Without windowing, positional coverage are outputted for the selected chromosomes.
///
/// ## Blacklisting
///
/// Positions in blacklisted regions are set to `f32::NaN` (and thus not included in sums or averages).
///
/// ## GC correction
///
/// Weight the contribution of each fragment by its length and GC content, using a precomputed
/// correction matrix (`cfdna gc-bias`). This reduces the global GC bias in the coverage,
/// which is a common technically-induced bias.
///
/// The GC correction matrix should be calculated from the same BAM file, as the bias is sample-specific.
///
/// ## Temporary files
///
/// We write temporary files to a `<output-dir>/tmp.<output-prefix>.<random>` directory to reduce memory.
/// This directory is deleted at the end of the run. If the software is disrupted, the directory
/// may be left behind.
///
/// ## Always-on exclusion criteria
///
/// The following criteria always exclude a read:
///
/// The read or mate read is unmapped.
/// The read is mapped to a different `tid` than the mate.
/// The read is secondary, supplementary or duplicate.
/// The read failed quality check.
/// The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct FCoverageConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.per_position.bedgraph.zst`,
    ///   `<prefix>.per_position_per_window.tsv.zst`,
    ///   `<prefix>.avg.tsv.zst`, or
    ///   `<prefix>.total.tsv.zst`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value = "coverage", help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Decimals to round coverage to when writing `[integer]`
    ///
    /// **NOTE**: When floating point precision is not needed,
    /// all coverages are integers, we remove all decimal points!
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "2", value_parser = clap::value_parser!(u8).range(0..), help_heading="Core"))]
    pub decimals: u8,

    /// Output zero-coverage runs in positional coverage outputs `[flag]`
    ///
    /// By default, only covered positions are written to the output.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub keep_zero_runs: bool,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    /// What to return per window `[string]`
    ///
    /// Possible values:
    ///
    /// - `"average"`: Get the average coverage per window (default).
    ///
    /// - `"total"`: Get the total coverage per window.
    ///
    /// - `"unique-positions"`: Get the positional coverage for the included windows only (`--by-bed` *only*).
    ///   Overlapping windows are merged to avoid duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"indexed-positions"`: Get the positional coverage for the included windows only (`--by-bed` *only*).
    ///   Adds the original window index as an output column and keeps duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// **NOTE**: Ignored when no windows are specified.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "average",
            value_parser,
            ignore_case = true,
            help_heading = "Core"
        )
    )]
    pub per_window: CoverageWindowAction,

    /// Ignore inter-mate gap `[flag]`
    ///
    /// Disable counting of the gap between reads (i.e., `[forward.end, reverse.start)`)
    /// when the two reads do not overlap.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub ignore_gap: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    ///
    /// Not recommended, as we already select only inward-directed read pairs within fragment length bounds.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    // TODO: Consider whether blacklist is "filtering" in tools like this?
    /// Optional BED file(s) with blacklisted regions `[path]`
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

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
}

impl FCoverageConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            output_prefix: "coverage".into(),
            decimals: 2,
            keep_zero_runs: false,
            tile_size: 20_000_000,
            per_window: CoverageWindowAction::Average,
            ignore_gap: false,
            windows: WindowsArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs {
                min_fragment_length: 20,
                max_fragment_length: 1000,
            },
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            gc: ApplyGCArgs { gc_file: None },
            ref_2bit: None,
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_decimals(&mut self, decimals: u8) {
        self.decimals = decimals;
    }

    pub fn set_scale_genome(&mut self, scale_genome: ScaleGenomeArgs) {
        self.scale_genome = scale_genome;
    }

    pub fn set_keep_zero_runs(&mut self, keep: bool) {
        self.keep_zero_runs = keep;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    pub fn set_per_window(&mut self, action: CoverageWindowAction) {
        self.per_window = action;
    }

    pub fn set_ignore_gap(&mut self, ignore: bool) {
        self.ignore_gap = ignore;
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
        self.windows = windows;
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}
