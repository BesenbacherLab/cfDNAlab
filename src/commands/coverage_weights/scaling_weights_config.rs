use std::path::PathBuf;

use crate::commands::cli_common::*;
use anyhow::bail;

/// Shared arguments for scaling-weight commands.
///
/// These settings control how raw fragment mass is counted before the
/// command-specific smoothing and scaling-factor postprocessing is applied.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct ScalingWeightsArgs {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'x',
            default_value_t = String::new(),
            hide_default_value = true,
            value_parser = crate::commands::cli_common::parse_output_prefix,
            help_heading = "Core"
        )
    )]
    pub output_prefix: String,

    /// Size (bp) of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "10000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    /// Size (bp) of the stride bins [integer]
    ///
    /// **NOTE**: `--bin-size` must be divisible by `stride`. I.e., `--bin-size % stride` == 0`.
    ///
    /// A multiplicative scaling factor is calculated per stride bin based on all its
    /// overlapping large-scale bins (`--bin-size`).
    ///
    /// Setting smaller values leads to a higher precision in the downstream normalization
    /// but also requires saving a larger TSV (one line per stride-bin)
    /// and takes longer to compute.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "500000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Core"))]
    pub stride: u32,

    /// Size (bp) of the large genomic bins used to build the triangularly weighted average [integer]
    ///
    /// Each stride bin is smoothed based on all large genomic bins that overlap it.
    /// Larger values lead to more smoothing across the genome as each
    /// stride bin is overlapped by more large bins from a broader region.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "5000000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Core"))]
    pub bin_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads [flag]
    ///
    /// This is **NOT** recommended by default, as it trims the tails of the length distribution.
    ///
    /// Note, that we only keep inward-directed fragments within the specified length range, so
    /// there's no real need for proper-pair filtering.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions [path]
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub gc: ApplyGCArgs,

    /// Optional 2bit reference genome file [path]
    ///
    /// NOTE: Required for `--gc-file`, otherwise ignored.
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

impl ScalingWeightsArgs {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            logging: LoggingArgs::default(),
            output_prefix: String::new(),
            bin_size: 5_000_000,
            stride: 500_000,
            tile_size: 10_000_000,
            chromosomes,
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

    pub fn set_output_prefix(&mut self, output_prefix: String) {
        self.output_prefix = output_prefix;
    }

    pub fn set_bin_size(&mut self, bin_size: u32) {
        self.bin_size = bin_size;
    }

    pub fn set_stride(&mut self, stride: u32) {
        self.stride = stride;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
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

    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }

    pub fn check_bin_sizes(&self) -> anyhow::Result<()> {
        let stride = self.stride;
        let bin_size = self.bin_size;

        if stride > bin_size {
            bail!(
                "stride ({}) cannot be higher than bin_size ({})",
                stride,
                bin_size,
            );
        }
        if !bin_size.is_multiple_of(stride) {
            bail!(
                "bin_size ({}) must be divisible by stride ({})",
                bin_size,
                stride
            );
        }

        Ok(())
    }
}
