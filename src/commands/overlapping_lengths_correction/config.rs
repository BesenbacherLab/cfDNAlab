//! Command-line and programmatic configuration for overlap-length model fitting.

use std::path::PathBuf;

use crate::commands::cli_common::{
    ChromosomeArgs, FragmentLengthArgs, IOCArgs, LoggingArgs, TempDirArgs, UnpairedArgs,
};
use crate::{ToCliCommand, cli_command::helpers::*};

/// Default number of core bases processed by each parallel tile.
pub const DEFAULT_TILE_SIZE: u32 = 5_000_000;
/// Lower fragment length limit used by the LIONHEART-derived model.
pub const DEFAULT_MIN_FRAGMENT_LENGTH: u32 = 100;
/// Inclusive upper fragment length limit used by the LIONHEART-derived model.
pub const DEFAULT_MAX_FRAGMENT_LENGTH: u32 = 220;
/// Default width of average overlapping fragment length bins.
pub const DEFAULT_LENGTH_BIN_SIZE: u32 = 3;
/// Default minimum mapping quality used by cfDNAlab cfDNA commands.
pub const DEFAULT_MIN_MAPQ: u8 = 30;

/// Fit LIONHEART's average overlapping fragment length normalization model.
///
/// The main purpose of this command is to let `fcoverage` produce coverage normalized for average
/// overlapping fragment length, as needed when reproducing LIONHEART scores. LIONHEART scores are
/// correlations between normalized cfDNA coverage and accessible-chromatin annotations. This
/// command does not calculate those correlations or write normalized coverage. It writes a model
/// package that can be applied with `fcoverage --overlap-length-file`.
///
/// At every eligible covered reference position, the command calculates the average fragment
/// length among the fragments covering that position. For paired fragments, cfDNAlab defines the
/// fragment length from `forward.pos` to `reverse.reference_end`. Positions with similar average
/// overlapping fragment lengths are grouped together, and the command estimates how their mean
/// coverage changes across the configured fragment length range.
///
/// The underlying hypothesis is that coverage can vary systematically with the lengths of the
/// fragments overlapping a position. Differences in sample fragmentation, library preparation,
/// or sequencing can therefore create coverage patterns that are unrelated to the chromatin
/// accessibility signal being measured. Without homogenizing this relationship, these patterns
/// can confound coverage-based LIONHEART scores and other analyses that use local fragment
/// coverage.
///
/// Positions covered by more fragments have a more stable average fragment length than positions
/// with low coverage. LIONHEART models this by letting the expected variation decrease with the
/// square root of raw fragment depth and combining the depth-specific distributions in a skewed
/// Student-t mixture. The first fit captures the smooth relationship between coverage and average
/// overlapping fragment length. Small irregularities around that curve and its asymmetry vary
/// between samples, so LIONHEART treats them as unwanted variation when samples are compared and
/// divides them out. A second fit captures the remaining broad shape. Its spread is retained, but
/// the distribution is made symmetric with a mean of 166 bp. This gives samples a common reference
/// without making their coverage profiles flat. LIONHEART uses 166 bp to make samples comparable.
/// It does not imply that another sample mean is biologically wrong. The resulting three
/// normalization factors are stored as a single multiplicative lookup weight for `fcoverage`.
///
/// This method was developed and evaluated primarily for fragment lengths from 100 to 220 bp.
/// Using a substantially different range is experimental. Like LIONHEART, this command always
/// fits the relationship from uncorrected coverage. GC correction and genomic scaling remain
/// independent `fcoverage` operations that can be combined with the fitted model during
/// application.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, PartialEq)]
pub struct OverlappingLengthsCorrectionConfig {
    /// Input BAM, output directory, and worker-thread settings.
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    /// Optional location for temporary per-run files.
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub temp: TempDirArgs,

    /// Whether each accepted read should be treated as a complete fragment.
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for the normalization model package `[string]`
    ///
    /// The output is `<prefix>.overlap_length_model.zarr`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Size of tiles to process in parallel `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_TILE_SIZE, value_parser = clap::value_parser!(u32).range(1000000..), help_heading = "Core")
    )]
    pub tile_size: u32,

    /// Width of average overlapping fragment length bins `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_LENGTH_BIN_SIZE, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Binning")
    )]
    pub length_bin_size: u32,

    /// Ignore the inter-mate gap when defining covered positions `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub ignore_gap: bool,

    /// Chromosomes included in the genome-wide sufficient statistics.
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Minimum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_MIN_FRAGMENT_LENGTH, value_parser = clap::value_parser!(u32).range(10..), help_heading = "Filtering")
    )]
    pub min_fragment_length: u32,

    /// Maximum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_MAX_FRAGMENT_LENGTH, value_parser = clap::value_parser!(u32).range(10..), help_heading = "Filtering")
    )]
    pub max_fragment_length: u32,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value_t = DEFAULT_MIN_MAPQ, help_heading = "Filtering")
    )]
    pub min_mapq: u8,

    /// Only include properly paired reads `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering")
    )]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Logging destination and verbosity settings.
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl OverlappingLengthsCorrectionConfig {
    /// Construct a normalization model configuration with the command's scientific defaults.
    ///
    /// The default 100-220 bp range is intentional. It is the range for which the underlying
    /// LIONHEART method was developed and should not be widened without treating the result as
    /// experimental.
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            temp: TempDirArgs::default(),
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix: String::new(),
            tile_size: DEFAULT_TILE_SIZE,
            length_bin_size: DEFAULT_LENGTH_BIN_SIZE,
            ignore_gap: false,
            chromosomes,
            min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
            max_fragment_length: DEFAULT_MAX_FRAGMENT_LENGTH,
            min_mapq: DEFAULT_MIN_MAPQ,
            require_proper_pair: false,
            blacklist: None,
            logging: LoggingArgs::default(),
        }
    }

    /// Return the fragment length filter in the shared command representation.
    ///
    /// Keeping this conversion in one place ensures tiling and BAM filtering use the same
    /// inclusive limits as the command configuration.
    pub fn fragment_lengths(&self) -> FragmentLengthArgs {
        FragmentLengthArgs {
            min_fragment_length: self.min_fragment_length,
            max_fragment_length: self.max_fragment_length,
        }
    }
}

impl ToCliCommand for OverlappingLengthsCorrectionConfig {
    /// Reconstruct the equivalent `cfdna overlap-length-model` command-line arguments.
    fn to_cli_args(&self) -> crate::Result<Vec<std::ffi::OsString>> {
        let mut args = command_args("overlap-length-model");
        push_ioc(&mut args, &self.ioc);
        push_temp_dir(&mut args, &self.temp);
        push_unpaired(&mut args, &self.unpaired);
        push_output_prefix(&mut args, &self.output_prefix);
        push_value(&mut args, "--tile-size", self.tile_size);
        push_value(&mut args, "--length-bin-size", self.length_bin_size);
        push_bool(&mut args, "--ignore-gap", self.ignore_gap);
        push_chromosomes(&mut args, &self.chromosomes);
        push_value(&mut args, "--min-fragment-length", self.min_fragment_length);
        push_value(&mut args, "--max-fragment-length", self.max_fragment_length);
        push_value(&mut args, "--min-mapq", self.min_mapq);
        push_bool(&mut args, "--require-proper-pair", self.require_proper_pair);
        push_path_values(&mut args, "--blacklist", self.blacklist.as_deref());
        push_logging(&mut args, &self.logging);
        Ok(args)
    }
}
