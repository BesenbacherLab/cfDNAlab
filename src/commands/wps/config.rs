use crate::commands::cli_common::{ApplyGCArgs, ScaleGenomeArgs};
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs, UnpairedArgs, WindowsArgs};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use std::path::PathBuf;

/// Calculate positional windowed protection scores (WPS) across the genome.
///
/// **Experimental**: enable via `--features cmd_wps` during `cargo build/install`.
///
/// In paired-end mode, only fragments with both reads present are considered.
///
/// NOTE: To extract nucleosome peaks via WPS, see `cfdna wps-peaks` instead.
///
/// WPS: Number of fragments fully overlapping the window, minus the number of fragments ending strictly inside the window.
/// Fragments that both start and end at the exact window edges are considered fully overlapping.
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
///  - Get the average or total WPS per window.
///
/// Without windowing, positional WPS are outputted for the selected chromosomes.
///
/// ## Blacklisting
///
/// Positions where the `--window_size` window overlaps a (dilated) blacklisted region are set to `f32::NaN` (and thus not included in sums or averages).
///
/// **Dilation**: We want to avoid any WPS scores being biased by neighbouring blacklisted intervals,
/// which can have an unreasonably high number of overlapping fragments.
/// Hence, we increase all blacklist intervals by the maximum fragment length + half the `--window_size` on both sides.
///
/// ## Scaling
///
/// When `--scaling-factors` are provided, we scale the **final per-position WPS** by the factor assigned
/// to the center base of that position.
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
/// The read is secondary, supplementary or duplicate.
/// The read failed quality check.
///
/// **Paired-end input only**:
/// The read or mate read is unmapped.
/// The read is mapped to a different `tid` than the mate.
/// The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).
///
/// ## Examples
///
/// ```rust,ignore
///
/// // Extract WPS scores (these arguments are always specified, hence `...` below)
/// cfdna wps --bam <> --output-dir <> -n-threads <>
///
/// // Extract positional WPS in windows
/// cfdna wps-peaks ... --by-bed <> --per-window "unique-positions"
///
/// ```
///
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct WPSConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub shared_args: WPSSharedConfig,

    /// Output zero-WPS runs in positional outputs `[flag]`
    ///
    /// By default, only positions with non-zero values are written to the output.
    #[cfg_attr(
        feature = "cli",
        clap(long, requires = "save-wps", help_heading = "Core")
    )]
    pub keep_zero_runs: bool,

    // TODO: For WPS, perhaps coefficient of variation is the relevant metric for window aggregation? std/mean ish?
    /// What to return for WPS per window `[string]`
    ///
    /// Possible values:
    ///
    /// - `"unique-positions"`: Get the positional WPS for the included windows only (`--by-bed` *only*).
    ///   Overlapping windows are merged to avoid duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"indexed-positions"`: Get the positional WPS for the included windows only (`--by-bed` *only*).
    ///   Adds the original window index as an output column and keeps duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"average"`: Get the average WPS per window.
    ///
    /// - `"total"`: Get the total WPS per window.
    ///
    /// **NOTE**: Ignored when no windows are specified.
    /// Required when `--save-wps` and either `--by-bed` or `--by-size` are provided.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            requires = "save-wps",
            ignore_case = true,
            help_heading = "Core"
        )
    )]
    pub per_window: Option<CoverageWindowAction>,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct WPSSharedConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

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

    /// Decimals to round values to when writing `[integer]`
    ///
    /// **NOTE**: When floating point precision is not needed,
    /// all values are integers, why we remove all decimal points!
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "2", value_parser = clap::value_parser!(u8).range(0..), help_heading="Core"))]
    pub decimals: u8,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

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

impl WPSSharedConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs, output_prefix: &str) -> Self {
        Self {
            ioc,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix: output_prefix.into(),
            window_size: 120,
            decimals: 2,
            tile_size: 20_000_000,
            windows: WindowsArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            min_fragment_length: 120,
            max_fragment_length: 180,
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            gc: ApplyGCArgs {
                gc_file: None,
                gc_tag: None,
                drop_invalid_gc: false,
            },
            ref_2bit: None,
        }
    }

    pub fn set_output_prefix(&mut self, prefix: String) {
        self.output_prefix = prefix;
    }

    pub fn set_window_size(&mut self, window_size: u32) {
        self.window_size = window_size;
    }

    pub fn set_decimals(&mut self, decimals: u8) {
        self.decimals = decimals;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
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

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}

impl WPSConfig {
    pub fn new(
        ioc: IOCArgs,
        chromosomes: ChromosomeArgs,
        per_window: Option<CoverageWindowAction>,
    ) -> Self {
        Self {
            shared_args: WPSSharedConfig::new(ioc, chromosomes, "wps"),
            keep_zero_runs: false,
            per_window: per_window,
        }
    }

    pub fn set_output_prefix(&mut self, prefix: String) {
        self.shared_args.set_output_prefix(prefix);
    }

    pub fn set_window_size(&mut self, window_size: u32) {
        self.shared_args.set_window_size(window_size);
    }

    pub fn set_decimals(&mut self, decimals: u8) {
        self.shared_args.set_decimals(decimals);
    }

    pub fn set_keep_zero_runs(&mut self, keep: bool) {
        self.keep_zero_runs = keep;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.shared_args.set_tile_size(tile_size);
    }

    pub fn set_per_window(&mut self, action: Option<CoverageWindowAction>) {
        self.per_window = action;
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
        self.shared_args.set_windows(windows);
    }

    pub fn set_min_fragment_length(&mut self, min_fragment_length: u32) {
        self.shared_args
            .set_min_fragment_length(min_fragment_length);
    }

    pub fn set_max_fragment_length(&mut self, max_fragment_length: u32) {
        self.shared_args
            .set_max_fragment_length(max_fragment_length);
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.shared_args.set_min_mapq(min_mapq);
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.shared_args.set_require_proper_pair(require);
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.shared_args.set_scale_genome(scale);
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.shared_args.set_gc(gc);
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.shared_args.set_ref_2bit(ref_2bit);
    }
}
