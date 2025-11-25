use crate::{
    commands::cli_common::{ApplyGCArgs, ChromosomeArgs, IOCArgs, ScaleGenomeArgs},
    shared::blacklist::BlacklistStrategy,
};
use std::path::PathBuf;

/// Count positional fragment **midpoint** coverage in groups of genomic windows.
///
/// **Midpoints**: The center of the fragment span (`[end(reverse), start(forward)]`),
/// with ties (in even-sized windows) randomly assigned to either the left or right mid-position
/// to reduce rounding bias.
///
/// **Groups**: The coverage profiles are "collapsed" (summed per position) for all windows in a group.
/// E.g., groups can be transcription factors with windows being binding sites. We then
/// get the overall profile per transcription factor.
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
pub struct ProfileGroupsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.midpoint_profiles.npy`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value = "sites", help_heading = "Core")
    )]
    pub output_prefix: String,

    /// The grouped fixed-size intervals to count within `[path]`
    ///
    /// A BED file of genomic intervals and their respective group names.
    ///
    /// Must be sorted by the `chromosome` and `start` coordinates, and
    /// all intervals must have the same length.
    ///
    /// Sites with the same group name are collapsed to a single profile.
    ///
    /// Columns: `chromosome, start, end, group_name`.
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

    /// Edges of fragment length bins to count in `[path]`
    ///
    /// The last edge is *exclusive*.
    ///
    /// **NOTE**: Memory consumption increases linearly with the number of bins.
    ///
    /// Example: `--length-bins 20 80 150 220 500 1001` or `--length-bins {20..1001..10}` for `20 30 40 ... 1001`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32).range(1..),
            num_args = 2.., // At least two edges per occurrence
            default_values_t = [20_u32, 1001_u32],
            help_heading = "Core"
        )
    )]
    pub length_bins: Vec<u32>,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "63000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

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
    /// This is NOT recommended by default as it trims the tails of the length distribution.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// **NOTE**: It may be an advantage to instead remove intervals that overlap
    /// blacklisted regions from the BED file.
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

impl ProfileGroupsConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs, intervals: PathBuf) -> Self {
        Self {
            ioc,
            output_prefix: "sites".into(),
            intervals,
            length_bins: vec![20, 1001],
            tile_size: 63_000_000,
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
            gc: ApplyGCArgs { gc_file: None },
            ref_2bit: None,
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_length_bins(&mut self, bins: Vec<u32>) {
        self.length_bins = bins;
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
