use std::path::PathBuf;

use crate::commands::cli_common::*;
use anyhow::bail;

/// Shared arguments for scaling-weight commands.
///
/// These settings control how raw fragment support is counted before the
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
            help_heading = "Core"
        )
    )]
    pub output_prefix: String,

    /// Size (bp) of large genomic bins to calculate coverage in [integer]
    ///
    /// Larger values lead to a more smooth coverage across the genome.
    ///
    /// **NOTE**: The normalizing scaling factors are calculated per stride-sized overlap
    /// of these bins. Technically, we only count the support per stride-sized bin
    /// and then calculate the overlap with a triangular weighting scheme.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "5000000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub bin_size: u32,

    /// Size (bp) of stride [integer]
    ///
    /// **NOTE**: `--bin-size` must be divisible by `stride`. I.e., `--bin-size % stride` == 0`.
    ///
    /// A normalizing scaling factor is calculated per stride as the (inverse) weighted average support of the overlapping large-scale bins.
    ///
    /// Smaller values lead to a higher precision in the downstream normalization
    /// but also require saving a larger TSV in the end (one line per stride-bin)
    /// and take longer to compute.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "500000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub stride: u32,

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
