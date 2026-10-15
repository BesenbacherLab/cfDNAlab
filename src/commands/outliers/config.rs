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
pub const DEFAULT_PROTECTION_WINDOW_SIZE: u32 = 100_000;
pub const DEFAULT_MIN_PROTECTION_POSITIONS: u32 = 50_000;
pub const DEFAULT_MAX_OUTLIER_PERCENTAGE: f64 = 5.0;

/// Marginal coverage distribution fitted in both stages of detection.
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutlierModel {
    /// Allow positive-count overdispersion, including the zero-dispersion ZIP boundary.
    #[default]
    Zinb,
    /// Fit the zero-inflated Poisson model without positive-count overdispersion.
    Zip,
}

impl OutlierModel {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Zinb => "zinb",
            Self::Zip => "zip",
        }
    }
}

/// Action shared by all calls in a merged region exceeding the preservation cutoff.
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AmbiguousRegionAction {
    /// Omit the merged region from correction weights and actionable blacklists.
    #[default]
    Preserve,
    /// Apply the proposed correction despite the density flag.
    Correct,
    /// Stop without publishing outputs when a merged region exceeds the cutoff.
    Error,
}

impl AmbiguousRegionAction {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Preserve => "preserve",
            Self::Correct => "correct",
            Self::Error => "error",
        }
    }
}

/// Coverage level retained after a region is called as an outlier.
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutlierTarget {
    /// Retain the active-component mean estimated by the second fit, excluding structural-zero dilution.
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

/// Detect smaller sample-specific pileups remaining after ordinary blacklisting.
///
/// The first tiled BAM pass builds raw positional coverage histograms in fixed cores. Neighboring
/// cores supply an adaptive local fitting context. A two-stage zero-inflated negative-binomial
/// (ZINB) model supplies thresholds, then a second BAM pass calls exact high-coverage intervals.
/// Use `--model zip` to compare the zero-inflated Poisson model with the same preservation rules.
///
/// Only high coverage is called. Eligible zeros remain in fitting and in window denominators.
/// Input-blacklisted positions are excluded from fitting, calling, and the automatic probability
/// denominator. Fragments crossing masked bases still contribute at their other covered positions.
///
/// ## Model and correction
///
/// ZINB combines structural zeros with a negative-binomial component having mean `mu` and variance
/// `mu + dispersion * mu^2`. Dispersion zero is ZIP. Both models first fit the complete histogram,
/// calculate T1, exclude coverage >= T1, and refit using the likelihood conditioned on X < T1.
/// T2 is the smallest positive integer with inclusive underlying tail P(X >= T2) at or below the
/// chosen probability. Preservation never adds excluded observations back to the baseline fit.
///
/// The automatic probability is 1 / eligible positions. Its interpretation is a fitted expected
/// position count, conditional on the model assumptions. It is not a family-wise error guarantee
/// and does not incorporate fitted-parameter uncertainty or model misspecification. Shared fragment
/// spans make adjacent coverage dependent. Inspect the observed-versus-fitted diagnostics.
///
/// Complete local contexts must fit successfully. Only a chromosome unable to meet the requested
/// span or fragment support uses the global fallback. Numerical failures report their context.
///
/// `local-mean` targets the fitted active-component mean, not the unconditional mean diluted by
/// structural zeros. Both means are reported. Each region's weight is its summed target coverage
/// divided by summed observed coverage, capped at 1. In `fcoverage --outlier-weights`, any positive
/// overlap applies the full regional weight to the complete original fragment. Multiple overlaps
/// use the minimum weight. A zero weight excludes the fragment before correction and counting.
///
/// ## Preserving broad or clustered elevations
///
/// Start a physical window at each exact call's start. A call uses the maximum called-base score
/// across all evaluated forward windows it overlaps, including earlier calls' windows. The score
/// is 100 * called eligible bases / max(eligible bases, --min-protection-positions).
///
/// If a call has no sufficiently supported forward window, try a window ending at its exact end.
/// Use that backward window if it has enough eligible bases. Otherwise keep the forward score
/// with the denominator floor. Windows are clipped at chromosome ends and never extended across
/// blacklisted gaps to collect eligible bases. The audit records the selected window and support.
///
/// A score strictly above `--max-outlier-percentage` flags the call. Flank every candidate, then
/// merge overlapping and touching flanks. Under `preserve`, any flagged member protects the whole
/// merged region. Its exact calls are withheld from weights and actionable BED files, and retained
/// in the audit. Flanking therefore affects shared decisions as well as blacklist extent.
///
/// This heuristic targets smaller residual pileups. It does not identify biological causes.
/// Narrow biological elevations can still be corrected, and dense technical events can be
/// preserved. The preservation window measures extent and does not supply a correction baseline.
///
/// ## Diagnostics and files
///
/// With an optional `<prefix>.`, the command writes:
///
/// - `outliers.keep_weights.tsv`: actionable sparse weights consumed by fcoverage
///
/// - `outliers.exact.bed` and `outliers.flanked.bed`: actionable exact and merged flanked blacklists
///
/// - `outliers.audit.tsv`: every candidate, proposed weight, window support, and final decision
///
/// - `outliers.summary.tsv`: per-chromosome and total call counts, bases, and eligible-base percentages
///
/// - `outliers.histograms.tsv.zst`: raw per-core coverage histograms and unblacklisted-midpoint fragment support
///
/// - `outliers.models.tsv`: local parameters, thresholds, context/fallback use, and fit diagnostics
///
/// - `outliers.global_fit.tsv` and `outliers.global_fit.png`: global distribution and local-threshold diagnostics
///
/// Positive-count variance ratios compare original observed positive values, including the tail,
/// with the chosen second-fit model's positive-count variance. Near 1 means comparable spread.
/// Tail counts at T1/T2 provide additional diagnostics. These values do not change calling.
/// All-zero incomplete contexts using global fallback have unavailable positive diagnostics (NaN).
///
/// Genomic call percentages do not measure the signal removed by complete-fragment weighting.
/// Compare unweighted and weighted fcoverage with the same filters and blacklists to inspect that
/// effect. Current support is limited to fcoverage. Consistent outlier-weighted GC/scaling
/// estimation is a future integration, not a prerequisite for multiplying existing GC weights.
///
/// Text outputs have `# key=value` settings comments. Coordinates are 0-based and half-open.
/// Rows follow selected chromosome order, then position. See the guide for an inspection workflow.
///
/// ## Fragment spans
///
/// Paired fragments span `forward.pos` to `reverse.reference_end`. Unpaired fragments span
/// `read.pos` to `read.reference_end`. Coverage respects mapped reference segments, so reference
/// D/N gaps do not contribute. `--ignore-gap` controls only the paired inter-mate gap.
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

    /// Minimum span (bp) used to fit each local coverage model `[integer]`
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
    /// called when its model survival probability `P(coverage >= observed)` is at most this value.
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
    /// `local-mean` uses the active-component mean from the underlying local model estimated
    /// by the second fit. `threshold` retains the final discrete calling
    /// threshold. `zero` excludes every fragment overlapping a called region when the weights are
    /// applied downstream.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_enum, default_value_t = OutlierTarget::LocalMean, help_heading = "Outlier Model")
    )]
    pub target: OutlierTarget,

    /// Coverage model used for the initial fit and right-truncated refit `[string]`
    ///
    /// ZINB allows extra variation in positive coverage. ZIP is available for direct comparisons.
    #[cfg_attr(feature = "cli", clap(long, value_enum, default_value_t = OutlierModel::Zinb, help_heading = "Outlier Model"))]
    pub model: OutlierModel,

    /// Physical reference span of each call-start preservation window `[integer]`
    ///
    /// Smaller windows preserve more localized clusters. This window does not define the fit.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_PROTECTION_WINDOW_SIZE, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Preservation"))]
    pub protection_window_size: u32,

    /// Minimum eligible-base denominator and support needed to avoid a backward fallback `[integer]`
    ///
    /// Short or masked windows use this floor. Low support alone does not preserve a region.
    /// This value must not exceed `--protection-window-size`.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_MIN_PROTECTION_POSITIONS, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Preservation"))]
    pub min_protection_positions: u32,

    /// Maximum called-base score before a merged region is flagged `[number]`
    ///
    /// Scores strictly above this percentage are flagged. Lower values preserve more regions.
    /// A value of 100 disables density-based preservation.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_MAX_OUTLIER_PERCENTAGE, help_heading = "Preservation"))]
    pub max_outlier_percentage: f64,

    /// Action for merged regions with a member above the percentage cutoff `[string]`
    #[cfg_attr(feature = "cli", clap(long, value_enum, default_value_t = AmbiguousRegionAction::Preserve, help_heading = "Preservation"))]
    pub ambiguous_region_action: AmbiguousRegionAction,

    /// Flank (bp) added to the separate blacklist BED `[integer]`
    ///
    /// When omitted, the maximum accepted fragment length is used. The exact unflanked BED and
    /// weighted interval TSV retain exact coordinates. Merged regions share preservation decisions.
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
            model: OutlierModel::Zinb,
            protection_window_size: DEFAULT_PROTECTION_WINDOW_SIZE,
            min_protection_positions: DEFAULT_MIN_PROTECTION_POSITIONS,
            max_outlier_percentage: DEFAULT_MAX_OUTLIER_PERCENTAGE,
            ambiguous_region_action: AmbiguousRegionAction::Preserve,
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

    /// Select the coverage model used in both fitting stages.
    pub fn set_model(&mut self, model: OutlierModel) {
        self.model = model;
    }

    /// Set the physical span of preservation windows.
    pub fn set_protection_window_size(&mut self, size: u32) {
        self.protection_window_size = size;
    }

    /// Set minimum eligible support and the short-window denominator floor.
    pub fn set_min_protection_positions(&mut self, positions: u32) {
        self.min_protection_positions = positions;
    }

    /// Set the strict percentage cutoff for shared merged-region decisions.
    pub fn set_max_outlier_percentage(&mut self, percentage: f64) {
        self.max_outlier_percentage = percentage;
    }

    /// Set the action applied to flagged merged regions.
    pub fn set_ambiguous_region_action(&mut self, action: AmbiguousRegionAction) {
        self.ambiguous_region_action = action;
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
        anyhow::ensure!(
            self.protection_window_size > 0,
            "protection_window_size must be positive"
        );
        anyhow::ensure!(
            self.min_protection_positions > 0
                && self.min_protection_positions <= self.protection_window_size,
            "min_protection_positions must be in 1..=protection_window_size, got {} and {}",
            self.min_protection_positions,
            self.protection_window_size
        );
        anyhow::ensure!(
            self.max_outlier_percentage.is_finite()
                && (0.0..=100.0).contains(&self.max_outlier_percentage),
            "max_outlier_percentage must be finite and in 0..=100, got {}",
            self.max_outlier_percentage
        );
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
        push_value(&mut args, "--model", self.model.as_str());
        push_value(
            &mut args,
            "--protection-window-size",
            self.protection_window_size,
        );
        push_value(
            &mut args,
            "--min-protection-positions",
            self.min_protection_positions,
        );
        push_value(
            &mut args,
            "--max-outlier-percentage",
            self.max_outlier_percentage,
        );
        push_value(
            &mut args,
            "--ambiguous-region-action",
            self.ambiguous_region_action.as_str(),
        );
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
