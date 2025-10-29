use std::path::PathBuf;

use crate::commands::cli_common::ScaleGenomeArgs;
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs, WindowsArgs};
use crate::commands::fcoverage::window_results::CoverageWindowAction;

/// Calculate positional windowed protections scores (WPS) across the genome.
///
/// Only paired-end fragments with both reads present are considered.
///
/// Tip: Use a minimum fragment length that matches the WPS `window_size`,
/// so no fragments have both ends in the WPS window.
///
/// ## Windowing (by-bed or by-size)
///
/// When specifying genomic windows via `--by-bed` or `--by-size`, one of the following outputs
/// is possible:
///
///  - Get the positional WPS for the included windows only (`--by-bed` *only*).
///    Excludes all positions that do not overlap a window from the output.
///    Choose between:
///     1) Indexed: Adds the original window index as an output column and keeps duplicate positions.
///     2) Unique: Overlapping windows are merged to avoid duplicate positions.
///
/// - Get the average or total WPS per window.
///
/// Without windowing, positional WPS are outputted for the selected chromosomes.
///
/// ## Blacklisting
///
/// Positions in blacklisted regions are set to `f32::NaN` (and thus not included in sums or averages).
/// Set `--nan-policy` to change how these positions are handled in the output (positional WPS outputs only).
///
/// ## Scaling
///
/// When `--scaling-factors` are provided, we scale the **final per-position WPS** by the factor assigned
/// to the centre base of that position.
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
pub struct WPSConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.wps.per_position.bedgraph.zst`,
    ///   `<prefix>.wps.per_position_per_window.tsv.zst`,
    ///   `<prefix>.wps.avg.tsv.zst`, or
    ///   `<prefix>.wps.total.tsv.zst`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value = "coverage", help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Window size `[integer]`
    ///
    /// Size of the window to calculate WPS scores from.
    /// Should be small enough that some fragments can fully overlap the window.
    #[cfg_attr(
        feature = "cli",
        clap(long, short='s', default_value = "120", value_parser = clap::value_parser!(u32).range(3..), help_heading="Core"))]
    pub window_size: u32,

    /// Decimals to round coverage to when writing `[integer]`
    ///
    /// **NOTE**: When floating point precision is not needed,
    /// all values are integers, why we remove all decimal points!
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "2", value_parser = clap::value_parser!(u8).range(0..), help_heading="Core"))]
    pub decimals: u8,

    /// Output zero-WPS runs in positional coverage outputs `[flag]`
    ///
    /// By default, only positions with non-zero values are written to the output.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub keep_zero_runs: bool,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    // TODO: Ensure a value is specified when windowing is enabled! indexed-positions would be the best default but that is --by-bed only?
    /// What to return per window `[string]`
    ///
    /// Possible values:
    ///
    /// - `"unique-positions"`: Get the positional coverage for the included windows only (`--by-bed` *only*).
    ///   Overlapping windows are merged to avoid duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"indexed-positions"`: Get the positional coverage for the included windows only (`--by-bed` *only*).
    ///   Adds the original window index as an output column and keeps duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"average"`: Get the average coverage per window.
    ///
    /// - `"total"`: Get the total coverage per window.
    ///
    /// **NOTE**: Ignored when no windows are specified.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, ignore_case = true, help_heading = "Core")
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

    /// Minimum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "120", value_parser = clap::value_parser!(u32).range(10..), help_heading="Filtering"))]
    pub min_fragment_length: u32,

    /// Maximum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "180", value_parser = clap::value_parser!(u32).range(10..), help_heading="Filtering"))]
    pub max_fragment_length: u32,

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
    // #[cfg_attr(feature = "cli", clap(flatten))]
    // gc: GCArgs,

    // #[cfg_attr(feature = "cli", clap(flatten))]
    // two_bit: TwoBitArgs,
}

impl WPSConfig {
    pub fn new(
        ioc: IOCArgs,
        chromosomes: ChromosomeArgs,
        per_window: CoverageWindowAction,
    ) -> Self {
        Self {
            ioc,
            output_prefix: "coverage".into(),
            window_size: 120,
            decimals: 2,
            keep_zero_runs: false,
            tile_size: 20_000_000,
            per_window: per_window,
            ignore_gap: false,
            windows: WindowsArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            min_fragment_length: 120,
            max_fragment_length: 180,
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_window_size(&mut self, window_size: u32) {
        self.window_size = window_size;
    }

    pub fn set_decimals(&mut self, decimals: u8) {
        self.decimals = decimals;
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

    pub fn set_min_fragment_length(&mut self, min_fragment_length: u32) {
        self.min_fragment_length = min_fragment_length;
    }

    pub fn set_max_fragment_length(&mut self, max_fragment_length: u32) {
        self.max_fragment_length = max_fragment_length;
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.scale_genome = scale;
    }
}
