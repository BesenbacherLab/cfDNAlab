use std::path::PathBuf;

use crate::commands::cli_common::{
    ChromosomeArgs, FragmentLengthArgs, IOCArgs, LoggingArgs, TempDirArgs, UnpairedArgs,
};
use crate::{ToCliCommand, cli_command::helpers::*};
use anyhow::{Result, bail};

pub const DEFAULT_STRIDE: u32 = 500_000;
pub const DEFAULT_BIN_SIZE: u32 = 5_000_000;
pub const DEFAULT_TILE_SIZE: u32 = 10_000_000;
pub const DEFAULT_MIN_CONTEXT_FRAGMENTS: u64 = 1_000;
pub const DEFAULT_TAIL_PROBABILITY_MULTIPLIER: f64 = 1.0;

/// Coverage level retained after a region is called as an outlier.
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutlierTarget {
    /// Retain `(1 - zero_inflation) * lambda` from the underlying ZIP estimated by the second fit.
    #[default]
    LocalMean,
    /// Retain coverage equal to the final discrete calling threshold.
    Threshold,
    /// Give every fragment overlapping the region a keep weight of zero.
    Zero,
}

impl OutlierTarget {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::LocalMean => "local-mean",
            Self::Threshold => "threshold",
            Self::Zero => "zero",
        }
    }
}

#[cfg(feature = "cli")]
fn parse_probability(value: &str) -> std::result::Result<f64, String> {
    let probability = value
        .parse::<f64>()
        .map_err(|_| format!("'{value}' is not a valid probability"))?;
    if !probability.is_finite() || probability <= 0.0 || probability > 0.5 {
        return Err("probability must be finite and in the interval (0, 0.5]".to_string());
    }
    Ok(probability)
}

#[cfg(feature = "cli")]
fn parse_positive_f64(value: &str) -> std::result::Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("'{value}' is not a valid number"))?;
    if !parsed.is_finite() || parsed <= 0.0 {
        return Err("value must be finite and greater than zero".to_string());
    }
    Ok(parsed)
}

/// Detect sample-specific regions with extreme raw fragment coverage.
///
/// The command counts uncorrected positional fragment coverage twice. The first tiled BAM pass
/// builds coverage histograms in `--stride`-sized cores. Neighboring histograms are combined until
/// each core has at least the requested `--bin-size` span and fragment support, then a two-stage
/// zero-inflated Poisson model supplies a local upper-tail threshold. The second pass calls the
/// exact outlier regions and calculates fragment keep weights.
///
/// Only unusually high coverage is called. Zeros and low coverage contribute to fitting but are
/// never themselves outliers. The default regional target is the expected local coverage from the
/// second ZIP fit, whose likelihood is conditioned on coverage below the first threshold.
///
/// ## Statistical model and limits
///
/// A ZIP variable is a mixture of a structural zero with probability `zero_inflation` and a
/// `Poisson(lambda)` value otherwise. Its expected coverage is
/// `(1 - zero_inflation) * lambda`. At positive coverage, zero inflation only multiplies the
/// Poisson probabilities by a constant. It cannot correct a positive-count distribution whose
/// shape or upper tail is not Poisson. Inspect the fitted-distribution outputs before treating ZIP
/// thresholds as calibrated for a new assay or sequencing depth.
///
/// The first threshold `T1` is fitted from the complete histogram. Coverage `>= T1` is excluded
/// from the second fit. That fit uses a right-truncated likelihood conditioned on coverage below
/// `T1`, but estimates the underlying untruncated ZIP. The final calling threshold `T2` is
/// calculated from the inclusive survival probability `P(X >= T2)` of that underlying model.
/// A core uses the global model only when its complete chromosome cannot supply the configured
/// context span or fragment support. If a complete local context cannot be fitted, the command
/// reports the error instead of silently substituting the global model.
///
/// The automatic probability `1 / number of eligible positions` limits the fitted sum of
/// position-wise tail probabilities to at most a single position. This is a plug-in expected-count
/// interpretation. It is not a family-wise error guarantee and does not account for uncertainty in
/// fitted parameters, model misspecification, or spatial dependence between positions covered by
/// the same fragments.
///
/// Outputs include sparse keep weights, exact and flanked BED files, per-core histograms with
/// fragment support, local model parameters, complete global fitted-distribution data, and a
/// diagnostic plot of both the global reference fit and the per-core models used during calling.
/// The keep-weight TSV is currently a standalone output. Other commands do not consume it yet.
///
/// ## Choosing an output
///
/// The two BED files make the calls usable immediately by commands and external tools that accept
/// BED blacklists. Use the exact BED when only the positions that crossed the outlier threshold
/// should be masked. Use the flanked BED for a more conservative blacklist that also masks nearby
/// positions potentially affected by fragments associated with the pileup. Blacklist behavior is
/// still determined by the command receiving the BED. The flank is an operational exclusion
/// buffer, not a claim that the added positions independently crossed the outlier threshold.
///
/// The keep-weight TSV is for the planned fractional fragment correction, where overlapping
/// fragments retain part of their contribution instead of being removed through blacklisting.
/// The histogram, model, global-fit, and plot outputs explain why regions were called and are for
/// model assessment rather than use as correction inputs.
///
/// ## Interpreting model diagnostics
///
/// `lambda` is the mean of the Poisson component, while `underlying_mean` is the expected coverage
/// after also accounting for structural zeros. `T1` excludes the initial extreme tail from the
/// second fit, and `T2` is the threshold actually used for calling.
///
/// The positive-coverage variance columns measure spread among positions whose coverage is greater
/// than zero. For an ordinary Poisson, variance equals `lambda`. After conditioning on positive
/// coverage, the expected variance is
/// `positive_mean * (1 + lambda - positive_mean)`. The reported ratio divides the variance of the
/// complete original positive coverage, including its extreme tail, by this expectation from the
/// underlying second-fit ZIP. A value near `1` is Poisson-like, a value above `1` is more variable,
/// and a value below `1` is less variable. It is a coarse diagnostic only and does not affect the
/// fitted model, thresholds, or calls. Observed positive metrics are `NaN` for an all-zero
/// incomplete context using the global fallback because it has no positive values to summarize.
///
/// The observed-versus-expected position counts at `T1` and `T2` are more direct checks of tail
/// calibration. Large discrepancies or consistently high variance ratios across representative
/// samples indicate that ZIP may not adequately describe positive coverage.
///
/// ## Output files
///
/// With an optional `<prefix>.`, the command writes:
///
/// - `outliers.keep_weights.tsv`: sparse fractional fragment weights for planned downstream use
///
/// - `outliers.exact.bed`: exact called positions for direct blacklisting or interval inspection
///
/// - `outliers.flanked.bed`: conservative blacklist with nearby positions added around each call
///
/// - `outliers.histograms.tsv.zst`: per-core raw data for auditing coverage and fragment support
///
/// - `outliers.models.tsv`: local fits, actual thresholds, fallback use, and fit diagnostics
///
/// - `outliers.global_fit.tsv`: machine-readable global observed-versus-fitted distribution
///
/// - `outliers.global_fit.png`: quick visual assessment of global fit and local thresholds
///
/// The public guide at https://cfdnalab.tools/docs/guides/handle_outliers_guide explains the
/// complete workflow and how to interpret every model metric.
///
/// Text outputs start with `# key=value` provenance comments. Coordinates are 0-based and
/// half-open. Rows follow the selected chromosome order and then genomic position.
///
/// ## Fragment span definition
///
/// **Paired-end**: `[forward.pos, reverse.reference_end)`, the reference span from the first
/// aligned position on the forward read to the last aligned position on the reverse read.
///
/// **Unpaired** where each read is a fragment: `[read.pos, read.reference_end)`.
///
/// Positional coverage follows the mapped reference segments of each read, so reference `D` and
/// `N` gaps do not contribute coverage. `--ignore-gap` controls only the inter-mate gap of a paired
/// fragment.
///
/// ## Blacklisting
///
/// Input-blacklisted positions do not influence the learned coverage expectations and are never
/// reported as outliers. Fragments that cross those positions still contribute at their other
/// covered positions.
///
/// ## Statistical references
///
/// - Lambert's ZIP definition and likelihood:
///   https://www.stat.cmu.edu/technometrics/90-00/vol-34-01/v3401001.pdf
///
/// - Stan's discrete truncation and endpoint definitions:
///   https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions
///
/// - NIST's inclusive discrete survival definition:
///   https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/intro.pdf
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, PartialEq)]
pub struct OutliersConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub temp: TempDirArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
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
    /// Chromosomes are processed in tiles of this size to limit peak memory use.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_TILE_SIZE, value_parser = clap::value_parser!(u32).range(1_000_000..), help_heading = "Core")
    )]
    pub tile_size: u32,

    /// Size (bp) of each fixed coverage-histogram core `[integer]`
    ///
    /// Every core receives a local model and final discrete coverage threshold. Smaller cores
    /// follow local baseline changes more closely but produce more model rows.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_STRIDE, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Core")
    )]
    pub stride: u32,

    /// Minimum span (bp) used to fit each local ZIP model `[integer]`
    ///
    /// Contexts expand by a complete `--stride` core on each available side. The realized span can
    /// therefore be larger than this value. With the defaults, centered contexts span 5.5 Mb.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_BIN_SIZE, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Core")
    )]
    pub bin_size: u32,

    /// Minimum accepted fragments required in a local fitting context `[integer]`
    ///
    /// Fragment support counts accepted fragments whose midpoint is not blacklisted, assigning
    /// each fragment to a single core. Contexts continue expanding after reaching `--bin-size`
    /// until this support is present. A core falls back to the global model only when its complete
    /// chromosome still lacks the configured span or support. A complete context that cannot be
    /// fitted is reported as an error.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_MIN_CONTEXT_FRAGMENTS, value_parser = clap::value_parser!(u64).range(1..), help_heading = "Core")
    )]
    pub min_context_fragments: u64,

    /// Manual upper-tail probability used to define outlier coverage `[number]`
    ///
    /// Without this option, the probability is `1 / number of eligible positions`. A position is
    /// called when its ZIP survival probability `P(coverage >= observed)` is at most this value.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser = parse_probability, help_heading = "Outlier Model")
    )]
    pub tail_probability: Option<f64>,

    /// Multiply the automatic upper-tail probability `[number]`
    ///
    /// Values above `1.0` lower the coverage threshold and call more positions. This multiplier is
    /// applied only to the automatic `1 / eligible positions` probability.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_TAIL_PROBABILITY_MULTIPLIER, value_parser = parse_positive_f64, conflicts_with = "tail_probability", help_heading = "Outlier Model")
    )]
    pub tail_probability_multiplier: f64,

    /// Coverage target retained inside called outlier regions `[string]`
    ///
    /// `local-mean` uses the expected coverage `(1 - zero_inflation) * lambda` from the underlying
    /// local ZIP estimated by the second fit. `threshold` retains the final discrete calling
    /// threshold. `zero` excludes every fragment overlapping a called region when the weights are
    /// applied downstream.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_enum, default_value_t = OutlierTarget::LocalMean, help_heading = "Outlier Model")
    )]
    pub target: OutlierTarget,

    /// Flank (bp) added to the separate blacklist BED `[integer]`
    ///
    /// When omitted, the maximum accepted fragment length is used. The exact unflanked BED and
    /// weighted interval TSV are unaffected.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser = clap::value_parser!(u32), help_heading = "Outlier Model")
    )]
    pub blacklist_flank: Option<u32>,

    /// Ignore the inter-mate gap when calculating raw coverage `[flag]`
    ///
    /// Cannot be used with `--reads-are-fragments`.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub ignore_gap: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading = "Filtering")
    )]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    ///
    /// This is not recommended because inward-directed fragments are already restricted by the
    /// fragment length bounds.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with positions excluded from fitting and calling `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering")
    )]
    pub blacklist: Option<Vec<PathBuf>>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl OutliersConfig {
    /// Create an outlier configuration with the documented command defaults.
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            temp: TempDirArgs::default(),
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix: String::new(),
            tile_size: DEFAULT_TILE_SIZE,
            stride: DEFAULT_STRIDE,
            bin_size: DEFAULT_BIN_SIZE,
            min_context_fragments: DEFAULT_MIN_CONTEXT_FRAGMENTS,
            tail_probability: None,
            tail_probability_multiplier: DEFAULT_TAIL_PROBABILITY_MULTIPLIER,
            target: OutlierTarget::LocalMean,
            blacklist_flank: None,
            ignore_gap: false,
            chromosomes,
            fragment_lengths: FragmentLengthArgs::default(),
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            logging: LoggingArgs::default(),
        }
    }

    /// Set the optional filename prefix placed before every output suffix.
    pub fn set_output_prefix<S: Into<String>>(&mut self, output_prefix: S) {
        self.output_prefix = output_prefix.into();
    }

    /// Set the directory used for temporary command files.
    pub fn set_temp_dir(&mut self, temp_dir: Option<PathBuf>) {
        self.temp.temp_dir = temp_dir;
    }

    /// Set the genomic tile size used to parallelize both BAM passes.
    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    /// Set the fixed core size used to assign local models and thresholds.
    pub fn set_stride(&mut self, stride: u32) {
        self.stride = stride;
    }

    /// Set the minimum genomic span of each adaptive local fitting context.
    pub fn set_bin_size(&mut self, bin_size: u32) {
        self.bin_size = bin_size;
    }

    /// Set the minimum midpoint-assigned fragment support for a local fitting context.
    pub fn set_min_context_fragments(&mut self, min_context_fragments: u64) {
        self.min_context_fragments = min_context_fragments;
    }

    /// Set a manual inclusive upper-tail probability, or `None` to use the automatic probability.
    pub fn set_tail_probability(&mut self, tail_probability: Option<f64>) {
        self.tail_probability = tail_probability;
    }

    /// Set the multiplier applied to the automatic `1 / eligible positions` probability.
    pub fn set_tail_probability_multiplier(&mut self, multiplier: f64) {
        self.tail_probability_multiplier = multiplier;
    }

    /// Set the coverage mass retained inside called regions.
    pub fn set_target(&mut self, target: OutlierTarget) {
        self.target = target;
    }

    /// Set the flank for the flanked BED, or use the maximum accepted fragment length by default.
    pub fn set_blacklist_flank(&mut self, blacklist_flank: Option<u32>) {
        self.blacklist_flank = blacklist_flank;
    }

    /// Control whether positional coverage excludes the inter-mate gap of paired fragments.
    pub fn set_ignore_gap(&mut self, ignore_gap: bool) {
        self.ignore_gap = ignore_gap;
    }

    /// Mutably access the accepted fragment length range.
    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    /// Set the minimum mapping quality accepted by the fragment iterator.
    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    /// Control whether paired-end reads must carry the proper-pair flag.
    pub fn set_require_proper_pair(&mut self, require_proper_pair: bool) {
        self.require_proper_pair = require_proper_pair;
    }

    /// Set BED files whose positions are excluded from fitting and calling.
    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

    pub(crate) fn validate(&self) -> Result<()> {
        self.fragment_lengths.validate()?;
        if self.tile_size < 1_000_000 {
            bail!("tile_size must be at least 1000000, got {}", self.tile_size);
        }
        if self.stride == 0 {
            bail!("stride must be greater than zero");
        }
        if self.bin_size == 0 {
            bail!("bin_size must be greater than zero");
        }
        if self.stride > self.bin_size {
            bail!(
                "stride ({}) cannot be greater than bin_size ({})",
                self.stride,
                self.bin_size
            );
        }
        if self.unpaired.reads_are_fragments && self.require_proper_pair {
            bail!("--require-proper-pair cannot be used with --reads-are-fragments");
        }
        if self.unpaired.reads_are_fragments && self.ignore_gap {
            bail!("--ignore-gap cannot be used with --reads-are-fragments");
        }
        if self.min_context_fragments == 0 {
            bail!("min_context_fragments must be greater than zero");
        }
        if !self.tail_probability_multiplier.is_finite() || self.tail_probability_multiplier <= 0.0
        {
            bail!(
                "tail_probability_multiplier must be finite and greater than zero, got {}",
                self.tail_probability_multiplier
            );
        }
        if let Some(probability) = self.tail_probability {
            if self.tail_probability_multiplier != DEFAULT_TAIL_PROBABILITY_MULTIPLIER {
                bail!(
                    "tail_probability and a non-default tail_probability_multiplier cannot be used together"
                );
            }
            if !probability.is_finite() || probability <= 0.0 || probability > 0.5 {
                bail!(
                    "tail_probability must be finite and in (0, 0.5], got {}",
                    probability
                );
            }
        }
        Ok(())
    }
}

impl ToCliCommand for OutliersConfig {
    fn to_cli_args(&self) -> crate::Result<Vec<std::ffi::OsString>> {
        let mut args = command_args("outliers");
        push_ioc(&mut args, &self.ioc);
        push_temp_dir(&mut args, &self.temp);
        push_unpaired(&mut args, &self.unpaired);
        push_output_prefix(&mut args, &self.output_prefix);
        push_value(&mut args, "--tile-size", self.tile_size);
        push_value(&mut args, "--stride", self.stride);
        push_value(&mut args, "--bin-size", self.bin_size);
        push_value(
            &mut args,
            "--min-context-fragments",
            self.min_context_fragments,
        );
        if let Some(probability) = self.tail_probability {
            push_value(&mut args, "--tail-probability", probability);
        } else if self.tail_probability_multiplier != DEFAULT_TAIL_PROBABILITY_MULTIPLIER {
            push_value(
                &mut args,
                "--tail-probability-multiplier",
                self.tail_probability_multiplier,
            );
        }
        push_value(&mut args, "--target", self.target.as_str());
        if let Some(flank) = self.blacklist_flank {
            push_value(&mut args, "--blacklist-flank", flank);
        }
        push_chromosomes(&mut args, &self.chromosomes);
        push_fragment_lengths(&mut args, &self.fragment_lengths);
        push_value(&mut args, "--min-mapq", self.min_mapq);
        push_bool(&mut args, "--require-proper-pair", self.require_proper_pair);
        push_path_values(&mut args, "--blacklist", self.blacklist.as_deref());
        push_logging(&mut args, &self.logging);
        push_bool(&mut args, "--ignore-gap", self.ignore_gap);
        Ok(args)
    }
}

#[cfg(test)]
mod tests {
    include!("config_tests.rs");
}
