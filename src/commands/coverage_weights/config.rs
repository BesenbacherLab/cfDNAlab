use crate::{commands::cli_common::*, shared::blacklist::BlacklistStrategy};
use anyhow::bail;
use std::path::PathBuf;

// TODO: Improve docstring - hard to understand

/// Extract fragment coverage in large genomic bins ("megabins") with a rolling
/// window and calculate normalizing scaling factors for smoothing the genome.
///
/// Outputs scaling factors per stride to allow other methods to apply the normalization (by weighting fragment counts).
///
/// The scaling factors are *inverted*, so normalization becomes multiplication.
/// Zero-valued coverages lead to zero-valued scaling factors. Non-zero factors have `mean == 1.0`.
///
/// ## Coverage
///
/// The full fragment span is counted without consideration of deletions and gaps.
/// This is fine for genome-scale normalization that reduces relative changes in coverage across the genome.
///
/// ## Fragment span
///
/// For **paired-end** sequencing, the span is defined as `[forward.pos, reverse.end)`.
/// For **unpaired** sequencing where each read is a fragment, the length is defined as `[read.pos, read.end)`.
///
/// ## Smoothing
///
/// Smoothing is performed as a triangular moving average, calculating
/// a weighted average of coverages from all bins overlapping a stride.  
///
/// ### Example
///
/// Assuming a bin-size of 6 and stride size of 2 (normally defaults to 5Mb and 0.5Mb respectively).
///
/// **Stride bins** (fixed along genome, each with an average coverage):
///
/// `[A] [B] [C] [D] [E] [F] [G] ...`
///
/// **Overlapping megabins** (`MB*`) (each covers 3 stride-bins). **`W_D`**, the number of overlapping megabins,
/// is the (unnormalized) weight of each stride-bin in the weighted-average coverage for stride-bin `D`:
///
/// ```text
///
/// MB1: [A][B][C]
///
/// MB2:    [B][C][D]
///
/// MB3:       [C][D][E]
///
/// MB4:          [D][E][F]
///
/// MB5:             [E][F][G]
///
/// W_D: [0][1][2][3][2][1][0]
///
/// ```
///
/// At chromosome edges, the weights are truncated (e.g., `W_D: [2][3][2][1][0]`).
///
/// The weights are normalized by their sum (after potential truncation at edges).
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
pub struct CoverageWeightsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.scaling_factors.tsv`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'x',
            default_value = "normalize_genome",
            help_heading = "Core"
        )
    )]
    pub output_prefix: String,

    /// Size (bp) of large genomic bins to calculate coverage in [integer]
    ///
    /// Larger values lead to a more smooth coverage across the genome.
    ///
    /// **NOTE**: The normalizing scaling factors are calculated per stride-sized overlap
    /// of these bins. Technically, we only count the coverage per stride-sized bin
    /// and then calculate the overlap with a triangular weighting scheme.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "5000000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub bin_size: u32,

    /// Size (bp) of stride [integer]
    ///
    /// **NOTE**: `--bin_size` must be divisible by `stride`. I.e., `bin_size % stride` == 0`.
    ///
    /// A normalizing scaling factor is calculated per stride as the (inverse) weighted average coverage of the overlapping large-scale bins.
    ///
    /// Smaller values lead to a higher precision in the downstream normalization
    /// but also require saving a larger BED file in the end (one line per stride-bin)
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

    /// Minimum size of blacklist intervals to load (bp) [integer]
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

    /// The fragment positions that should overlap blacklisted regions for it to be excluded [string]
    ///
    /// Possible values:
    ///     `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"` [string]
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
}

impl CoverageWeightsConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix: "normalize_genome".to_string(),
            bin_size: 5_000_000,
            stride: 500_000,
            chromosomes,
            fragment_lengths: FragmentLengthArgs::default(),
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
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

    pub fn check_bin_sizes(&self) -> anyhow::Result<()> {
        let stride = self.stride.clone();
        let bin_size = self.bin_size.clone();

        if stride > bin_size {
            bail!(
                "stride ({}) cannot be higher than bin_size ({})",
                stride,
                bin_size,
            );
        }
        if bin_size % stride != 0 {
            bail!(
                "bin_size ({}) must be divisible by stride ({})",
                bin_size,
                stride
            );
        }

        Ok(())
    }
}
