use std::path::PathBuf;

use crate::commands::cli_common::{ApplyGCArgs, ScaleGenomeArgs};
use crate::commands::cli_common::{
    ChromosomeArgs, FragmentLengthArgs, IOCArgs, LoggingArgs, UnpairedArgs, WindowsArgs,
};
use crate::commands::fcoverage::window_results::CoverageWindowAction;

/// Count positional **fragment** coverage across the genome.
///
/// In paired-end mode, only fragments with both reads present are considered.
/// By default, the entire fragment span is counted, except for
/// deletions and skipped regions that are not covered by the other read.
///
/// ## Fragment span definition
///
/// **Paired-end**: `[forward.pos, reverse.end)`.
///
/// **Unpaired** where each read is a fragment: `[read.pos, read.end)`.
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
/// ## Positional output and tiles
///
/// **Positional** outputs are written tile by tile to keep memory use low.
/// This means coverage segments can be split at genomic tile boundaries even when the
/// coverage value stays the same.
/// The covered positions and coverage values stay the same, but the bedGraph rows
/// may be shorter than they would be in a single-pass whole-chromosome run.
///
/// Reduced outputs like per-window `average` and `total` are merged across tiles,
/// so tile boundaries should not affect their final values.
///
/// ## Blacklisting
///
/// Positions in blacklisted regions are set to `f32::NaN` (and thus not included in sums or averages).
///
/// ## GC correction
///
/// Reduce the global GC bias (common technically-induced bias) in the coverage
/// by weighting the contribution of fragments. Two options:
///
/// `--gc-file`: Weight the contribution of each fragment by its length and GC content using a precomputed
/// correction matrix from `cfdna gc-bias`. The GC correction matrix should be calculated from the same BAM file,
/// as the bias is sample-specific.
///
/// `--gc-tag`: Weight the contribution of each fragment by a weight saved as an aux tag in the BAM reads.
/// Allows using external GC packages like `GCParagon` and `GCfix` (both use the tag "GC").
///
/// ## Temporary files
///
/// We write temporary files to a `<output-dir>/tmp.<output-prefix>.<random>` directory to reduce memory.
/// When no output prefix is given, the directory becomes `<output-dir>/tmp.<random>`.
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
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct FCoverageConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Divide the contribution of each fragment by the number of countable bases [flag]
    ///
    /// By default, we count each fragment as `1.0` in each covered position (before correction/scaling).
    /// That weights longer fragments higher than shorter fragments in the overall mass
    /// as they are counted in more positions. If we want each fragment
    /// to contribute the **same mass**, we can divide the per-position
    /// `1.0` weight by the number of countable positions.
    ///
    /// **Interpretation**: Per-base fragment support after normalizing each fragment
    /// to a total weight of `1.0` before correction/scaling.
    /// For `--per-window total` this approximates fragment counts
    /// (in sufficiently large windows).
    ///
    /// This flag is reflected in the output filenames.
    ///
    /// Blacklisted positions still count toward the normalization denominator
    /// to avoid large values around blacklisted regions (edge effects).
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub normalize_by_length: bool,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.fcoverage.per_position.bedgraph.zst`,
    ///   `<prefix>.fcoverage.per_position_per_window.tsv.zst`,
    ///   `<prefix>.fcoverage.avg.tsv.zst`, or
    ///   `<prefix>.fcoverage.total.tsv.zst`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, help_heading = "Core")
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
        clap(long, default_value = "10000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
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
    ///   **NOTE**: The output is first sorted by chromosome, tile index, and window start.
    ///   Then the coverage segments are sorted by start- and end coordinates.
    ///   Window indices may thus not be contiguous.
    ///   Depending on your needs, sort downstream.
    ///   
    ///   
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
    ///
    /// Cannot be used with `--reads-are-fragments`.
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

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl FCoverageConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            logging: LoggingArgs::default(),
            normalize_by_length: false,
            output_prefix: String::new(),
            decimals: 2,
            keep_zero_runs: false,
            tile_size: 10_000_000,
            per_window: CoverageWindowAction::Average,
            ignore_gap: false,
            windows: WindowsArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs::default(),
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            gc: ApplyGCArgs {
                gc_file: None,
                gc_tag: None,
                neutralize_invalid_gc: false,
            },
            ref_2bit: None,
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_unpaired(&mut self, unpaired: UnpairedArgs) {
        self.unpaired = unpaired;
    }

    pub fn set_normalize_by_length(&mut self, normalize_by_length: bool) {
        self.normalize_by_length = normalize_by_length;
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

    pub fn set_fragment_lengths(&mut self, fragment_lengths: FragmentLengthArgs) {
        self.fragment_lengths = fragment_lengths;
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}
