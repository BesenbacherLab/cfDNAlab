use std::{fmt, path::PathBuf, str::FromStr};

use crate::commands::cli_common::{
    ChromosomeArgs, FragmentLengthArgs, IOCArgs, LoggingArgs, TempDirArgs, UnpairedArgs,
};
use crate::{ToCliCommand, cli_command::helpers::*};
use anyhow::{Context, Result, bail, ensure};

pub const DEFAULT_STRIDE: u32 = 500_000;
pub const DEFAULT_BIN_SIZE: u32 = 5_000_000;
pub const DEFAULT_TILE_SIZE: u32 = 10_000_000;
pub const DEFAULT_MIN_CONTEXT_FRAGMENTS: u64 = 1_000;
pub const DEFAULT_TAIL_PROBABILITY_MULTIPLIER: f64 = 1.0;
pub const DEFAULT_PROTECTION_WINDOW_SIZE: u32 = 100_000;
pub const DEFAULT_MIN_PROTECTION_POSITIONS: u32 = 50_000;
pub const DEFAULT_MAX_OUTLIER_PERCENTAGE: f64 = 5.0;
pub const DEFAULT_PRESERVED_PLOT_CONTEXT_MULTIPLIER: f64 = 10.0;
pub const DEFAULT_PRESERVED_PLOT_MIN_SPAN: u32 = 100_000;

/// Statistical model used to estimate the expected range of fragment coverage.
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutlierModel {
    /// Zero-inflated negative binomial, allowing extra variation in coverage.
    #[default]
    Zinb,
    /// Zero-inflated Poisson, for comparison with a simpler coverage model.
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

/// How surrounding model cores contribute to a local coverage estimate.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ContextWeighting {
    /// Give every eligible base in the context equal influence.
    #[default]
    Uniform,
    /// Give nearby cores more influence, decreasing linearly with distance from the center core.
    Triangular,
    /// Blend normalized uniform and triangular histograms using this triangular fraction.
    TriangularBlend(f64),
}

impl ContextWeighting {
    /// Return the triangular contribution, rejecting invalid programmatic blend settings.
    pub(crate) fn triangular_fraction(self) -> Result<f64> {
        let fraction = match self {
            Self::Uniform => 0.0,
            Self::Triangular => 1.0,
            Self::TriangularBlend(fraction) => fraction,
        };
        ensure!(
            fraction.is_finite() && (0.0..=1.0).contains(&fraction),
            "context weighting blend fraction must be finite and in [0, 1], got {fraction}"
        );
        Ok(fraction)
    }
}

impl fmt::Display for ContextWeighting {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Uniform => formatter.write_str("uniform"),
            Self::Triangular => formatter.write_str("triangular"),
            Self::TriangularBlend(fraction) => write!(formatter, "triangular={fraction}"),
        }
    }
}

impl FromStr for ContextWeighting {
    type Err = anyhow::Error;

    fn from_str(argument: &str) -> Result<Self> {
        let weighting = match argument {
            "uniform" => Self::Uniform,
            "triangular" => Self::Triangular,
            _ => {
                let fraction = argument.strip_prefix("triangular=").context(
                    "use uniform, triangular, or triangular=<fraction> with a fraction in [0, 1]",
                )?;
                Self::TriangularBlend(
                    fraction
                        .parse()
                        .context("parsing context weighting blend fraction")?,
                )
            }
        };
        weighting.triangular_fraction()?;
        Ok(weighting)
    }
}

/// How to handle regions where too much nearby sequence is called as an outlier.
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AmbiguousRegionAction {
    /// Leave these regions uncorrected and omit them from the output weights and BED files.
    #[default]
    Preserve,
    /// Correct these regions in the same way as other outliers.
    Correct,
    /// Stop without publishing results if any region exceeds the percentage cutoff.
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

/// Coverage level used to calculate the correction for an outlier region.
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutlierTarget {
    /// Reduce coverage toward the estimated background mean.
    #[default]
    LocalMean,
    /// Reduce coverage toward the cutoff used to detect outliers.
    Threshold,
    /// Exclude every fragment overlapping the outlier region when applying the weights.
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

/// Parse the plot multiplier using the same bounds as programmatic validation.
#[cfg(feature = "cli")]
fn parse_context_multiplier(value: &str) -> std::result::Result<f64, String> {
    let parsed = parse_positive_f64(value)?;
    if parsed < 1.0 {
        return Err("plot context multiplier must be at least 1".to_string());
    }
    Ok(parsed)
}

/// Reject invalid preservation percentages before opening any input files.
#[cfg(feature = "cli")]
fn parse_percentage(value: &str) -> std::result::Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("'{value}' is not a valid percentage"))?;
    if !parsed.is_finite() || !(0.0..=100.0).contains(&parsed) {
        return Err("percentage must be finite and in 0..=100".to_string());
    }
    Ok(parsed)
}

/// Find sample-specific regions with unusually high fragment coverage and calculate correction weights.
///
/// Use this after excluding known problematic regions with your usual blacklist. It is intended
/// to reduce the influence of extreme pileups over short genomic intervals in individual samples.
///
/// The command writes weights for `fcoverage --outlier-weights` and BED files for excluding the
/// detected regions. It also reports how much of the genome was called and which regions were
/// left uncorrected. Inspect the results before adopting the correction in an analysis pipeline.
///
/// ## How outliers are found
///
/// Coverage is the number of fragments covering a genomic position, before GC correction or other
/// weighting. For each part of a chromosome, the command estimates the expected range of coverage
/// from a larger surrounding region. Positions at or above the resulting coverage cutoff are called as
/// outliers. Adjacent called positions are joined into regions.
///
/// The coverage model is fitted twice. After the first fit, the most extreme values are excluded
/// from the second fit so they have less influence on the final cutoff. Choose the statistical
/// model with `--model` and control how extreme coverage must be with `--tail-probability`.
///
/// Eligible positions are positions outside the supplied blacklists on the selected chromosomes.
/// They contribute to the estimates even when their coverage is zero. Only high coverage is called
/// as an outlier. Fragments crossing a blacklist can still contribute coverage outside it.
///
/// ## Applying the correction
///
/// Each selected outlier region gets a keep weight between 0 and 1. A weight of 0.1 retains 10%
/// of a fragment's contribution, while 0 excludes the fragment. By default, weights aim to reduce
/// the region's coverage toward the estimated background mean. See `--target` for alternatives.
///
/// `fcoverage` applies the full regional weight to every overlapping fragment, including its
/// contribution outside the outlier region. Even a 1 bp overlap applies the full weight. If a
/// fragment overlaps multiple regions, the smallest weight is used. Other fragments are unchanged.
/// Use the same fragment filters and blacklists when detecting outliers and applying the weights.
///
/// ## Broad coverage elevations
///
/// High coverage can also reflect biological changes, such as copy-number amplification. To help
/// avoid correcting broad elevations, the command checks how much of the surrounding sequence was
/// called as an outlier. Regions exceeding `--max-outlier-percentage` are left uncorrected unless
/// you choose a different `--ambiguous-region-action`. They remain listed in the audit file.
///
/// A window starts at each called region. Each region uses the highest outlier percentage from
/// any of these windows that overlaps it, including windows started by earlier regions. Window
/// size is set by `--protection-window-size`. Blacklisted bases are excluded from the percentage,
/// but zero-coverage positions are included. Windows below `--min-protection-positions` are ignored
/// whenever another overlapping window meets that minimum.
///
/// Calls are also extended on both sides by `--blacklist-flank`. Overlapping or touching extended
/// regions are joined. If any part is being left uncorrected, all calls in the joined region are
/// left uncorrected. The weights themselves still apply to the original, unextended calls.
///
/// This rule measures the extent of the coverage elevation. It cannot identify its biological or
/// technical cause. Small biological events may still be corrected, and broad technical pileups
/// may be left uncorrected. It does not change the surrounding coverage estimate.
///
/// ## Output files
///
/// An optional `--output-prefix` is followed by a dot before these filenames:
///
/// - `outliers.keep_weights.tsv`: weights to use with `fcoverage --outlier-weights`.
///
/// - `outliers.exact.bed` and `outliers.flanked.bed`: regions to exclude, without and with flanking.
///   Both omit the regions being left uncorrected.
///
/// - `outliers.audit.tsv`: every called region, its proposed weight, and the reason for correcting
///   it or leaving it uncorrected.
///
/// - `outliers.summary.tsv`: numbers of called regions and bases, split by whether they are
///   selected for correction or left unchanged. Percentages are relative to eligible positions.
///
/// - `outliers.preserved_regions.tsv`: preserved calls grouped for inspection, with their
///   approximate spans and plotting coordinates. Use `--plot-preserved-regions` to add PNG plots.
///
/// - `outliers.models.tsv`, `outliers.global_fit.tsv`, and `outliers.global_fit.png`: coverage
///   estimates, detection cutoffs, and a plot comparing observed coverage with the fitted models.
///
/// - `outliers.histograms.tsv.zst`: counts of positions at each coverage level, for each
///   `--stride` interval.
///
/// - `outliers.README.txt`: explanations of the files, columns, and preservation decisions.
///
/// Genomic intervals use 0-based start coordinates and exclude the end coordinate, as in BED files.
///
/// The percentage of called genomic bases differs from the fraction of coverage removed, because
/// weights follow complete fragments. Compare unweighted and weighted `fcoverage` results to
/// assess the effect on your analysis. The current fragment-weighting integration is in `fcoverage`.
///
/// ## Example
///
/// `cfdna outliers --bam sample.bam --output-dir outliers --blacklist blacklist.bed`
///
/// Apply the resulting weights with the same filters and blacklist:
///
/// `cfdna fcoverage --bam sample.bam --output-dir coverage --blacklist blacklist.bed --outlier-weights outliers/outliers.keep_weights.tsv`
///
/// ## Fragment span definition
///
/// **Paired-end**: `[forward.pos, reverse.reference_end)`, from the first aligned position on
/// the forward read to the last aligned position on the reverse read.
///
/// **Unpaired**, where each read is a fragment: `[read.pos, read.reference_end)`.
///
/// Deletions and skipped reference regions do not contribute coverage unless covered by the
/// other read. `--ignore-gap` additionally excludes the gap between paired reads.
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

    /// Spacing (bp) of local coverage estimates along each chromosome `[integer]`
    ///
    /// Split each chromosome into intervals of this size. Each interval uses a single coverage
    /// cutoff estimated from a larger surrounding region, whose eligible-base minimum is `--bin-size`.
    /// Smaller values update the cutoff more frequently but require more model fits.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_STRIDE, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Core")
    )]
    pub stride: u32,

    /// Minimum number of eligible bases used to estimate local coverage `[integer]`
    ///
    /// Start with the central interval defined by `--stride`, then add an interval on each side
    /// until their eligible bases reach this minimum. Include zero-coverage bases, but exclude
    /// blacklisted bases. The physical span can therefore be larger. Without blacklists, 500 kb
    /// intervals with a 5 Mb minimum give centered spans of 0.5, 1.5, 2.5, ..., 5.5 Mb.
    ///
    /// At chromosome ends, expansion continues on the available side. The region can grow further
    /// if needed to include `--min-context-fragments` fragments.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_BIN_SIZE, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Core")
    )]
    pub bin_size: u32,

    /// Weight surrounding coverage observations by distance from the central interval `[choice]`
    ///
    /// Triangular weighting gives the central `--stride` interval the most influence, decreasing
    /// linearly by whole intervals toward the context edges. Both fitting stages use these weights.
    /// This keeps the estimate local while borrowing information from a wider region.
    /// Uniform weighting gives every eligible base equal influence.
    ///
    /// Use `triangular=0.2` for an 80% uniform and 20% triangular blend. Normalize both
    /// histograms before blending, so masking does not change the requested proportions.
    /// Fractions range from 0 (uniform) to 1 (fully triangular).
    ///
    /// Coverage values and the eligible-base and fragment support minima are not weighted.
    /// At chromosome ends the triangle stays centered on the interval being estimated, using
    /// only the available side. Sample-wide fallback models are always fitted without weighting.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = ContextWeighting::Uniform, help_heading = "Core"))]
    pub context_weighting: ContextWeighting,

    /// Minimum number of fragments used to estimate local coverage `[integer]`
    ///
    /// Count fragments whose midpoint is inside the surrounding region and outside the blacklists.
    /// If there are too few, keep expanding after the eligible-base minimum is reached.
    ///
    /// If the whole chromosome cannot provide both the eligible-base minimum and fragment count, use
    /// the coverage model estimated from all selected chromosomes instead.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_MIN_CONTEXT_FRAGMENTS, value_parser = clap::value_parser!(u64).range(1..), help_heading = "Core")
    )]
    pub min_context_fragments: u64,

    /// Probability cutoff for calling unusually high coverage `[number]`
    ///
    /// Call a position when the fitted model gives this probability or less of seeing its
    /// coverage or higher. Use smaller values for a stricter cutoff that requires more extreme
    /// coverage, and larger values for a more permissive cutoff.
    ///
    /// When omitted, use `1 / number of eligible positions` across the selected chromosomes,
    /// including positions with zero coverage. The actual false-positive rate depends on how well
    /// the model describes the data. Values must be greater than 0 and at most 0.5.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser = parse_probability, help_heading = "Outlier Model")
    )]
    pub tail_probability: Option<f64>,

    /// Adjust the automatically chosen probability cutoff `[number]`
    ///
    /// Multiply the automatic cutoff by this value. Values above `1.0` make detection more
    /// permissive, while values below `1.0` make it stricter. Cannot be used with a manually
    /// supplied `--tail-probability`.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_TAIL_PROBABILITY_MULTIPLIER, value_parser = parse_positive_f64, conflicts_with = "tail_probability", help_heading = "Outlier Model")
    )]
    pub tail_probability_multiplier: f64,

    /// Coverage level used to calculate how much to downweight an outlier region `[string]`
    ///
    /// `local-mean` uses the estimated background mean. The model handles additional zero-coverage
    /// positions separately, so those extra zeros do not lower this target. Both the background
    /// mean and the mean including those zeros are reported in `outliers.models.tsv`.
    ///
    /// `threshold` uses the cutoff for detecting outliers as the target coverage. `zero`
    /// excludes every fragment overlapping the region when the weights are applied.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_enum, default_value_t = OutlierTarget::LocalMean, help_heading = "Outlier Model")
    )]
    pub target: OutlierTarget,

    /// Statistical model used to estimate the expected range of coverage `[string]`
    ///
    /// `zinb` is a zero-inflated negative binomial model. It allows coverage to vary more than
    /// a Poisson model. `zip` is a zero-inflated Poisson model, available for comparison.
    /// Both account for positions with zero coverage and use the same two fitting steps.
    #[cfg_attr(feature = "cli", clap(long, value_enum, default_value_t = OutlierModel::Zinb, help_heading = "Outlier Model"))]
    pub model: OutlierModel,

    /// Window size (bp) used to recognize broad or clustered coverage elevations `[integer]`
    ///
    /// Start a window of this size at each called region. Each region uses the highest percentage
    /// of outlier bases from these windows that overlap it, including earlier regions' windows.
    /// Use only windows meeting `--min-protection-positions` when any overlap the call.
    /// Smaller values focus this check on more localized elevations.
    ///
    /// Windows stop at chromosome ends. Blacklisted bases are excluded from the counts without
    /// extending the genomic span. This setting is separate from the coverage-estimation region
    /// controlled by `--bin-size`.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_PROTECTION_WINDOW_SIZE, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Preservation"))]
    pub protection_window_size: u32,

    /// Minimum number of eligible positions used for window percentages `[integer]`
    ///
    /// A window may have fewer usable bases near a chromosome end or a large blacklisted region.
    /// If any overlapping forward window meets this minimum, take the highest percentage from
    /// those windows and ignore windows with fewer bases.
    ///
    /// If none qualifies, also check a window ending at the call's end. Use it if it meets the
    /// minimum. Otherwise, take the highest percentage from the original windows, dividing each
    /// outlier-base count by at least this many bases.
    ///
    /// This prevents a very short window from producing a large percentage based on little
    /// sequence. Having too few bases alone does not keep a region uncorrected.
    /// Cannot exceed `--protection-window-size`.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_MIN_PROTECTION_POSITIONS, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Preservation"))]
    pub min_protection_positions: u32,

    /// Maximum percentage of outlier bases allowed in a surrounding window `[number]`
    ///
    /// Percentages strictly above this value trigger `--ambiguous-region-action`. With `preserve`,
    /// the affected regions are left uncorrected. Equality does not trigger the action.
    /// With `preserve`, lower values leave more regions uncorrected. A value of 100 disables this check.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_MAX_OUTLIER_PERCENTAGE, value_parser = parse_percentage, help_heading = "Preservation"))]
    pub max_outlier_percentage: f64,

    /// How to handle broad or clustered coverage elevations `[string]`
    ///
    /// Applies when a region's window exceeds `--max-outlier-percentage`. `preserve` leaves it
    /// uncorrected, `correct` treats it as an ordinary outlier, and `error` reports every failing
    /// region after the second BAM pass and stops without publishing normal results. A readable
    /// `outliers.error.txt` and an `outliers.error_audit.tsv` explain the failures.
    /// Regions joined after `--blacklist-flank` share the same action. Preserved calls remain in
    /// `outliers.audit.tsv` and are omitted from the weights and BED files.
    #[cfg_attr(feature = "cli", clap(long, value_enum, default_value_t = AmbiguousRegionAction::Preserve, help_heading = "Preservation"))]
    pub ambiguous_region_action: AmbiguousRegionAction,

    /// Plot coverage around each group of preserved calls `[flag]`
    ///
    /// Add a PNG with surrounding coverage and a closer view of each group listed in
    /// `outliers.preserved_regions.tsv`. Show raw coverage, the estimated background mean,
    /// the outlier cutoff, called positions, and blacklisted gaps.
    ///
    /// Plotting reads the BAM again only around these groups, using the detection filters.
    /// The table is always written. It groups preserved calls using overlapping preservation
    /// windows and reports approximate locations, not biological breakpoints. Grouping is for
    /// inspection and does not change correction.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Preserved Region Plots"))]
    pub plot_preserved_regions: bool,

    /// Multiply the grouped region's span to choose the total plot width `[number]`
    ///
    /// Center this wider view on the grouped calls. Include the preservation windows and at least
    /// `--preserved-plot-min-span` bases, and stop at chromosome bounds. Values must be at least 1.
    /// Plot coordinates are recorded in the table even when plots are not requested.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_PRESERVED_PLOT_CONTEXT_MULTIPLIER, value_parser = parse_context_multiplier, help_heading = "Preserved Region Plots"))]
    pub preserved_plot_context_multiplier: f64,

    /// Minimum total genomic span (bp) shown around a preserved group `[integer]`
    ///
    /// Give short calls enough surrounding coverage to inspect. This includes blacklisted gaps
    /// and is separate from the minimum eligible-position count used for preservation decisions.
    #[cfg_attr(feature = "cli", clap(long, default_value_t = DEFAULT_PRESERVED_PLOT_MIN_SPAN, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Preserved Region Plots"))]
    pub preserved_plot_min_span: u32,

    /// Extend outlier regions by this many bases on each side in the flanked BED `[integer]`
    ///
    /// When omitted, use the maximum accepted fragment length. Overlapping or touching extended
    /// regions are joined. If any call in a joined region is being left uncorrected, the other
    /// calls in that region are also left uncorrected.
    ///
    /// Weights and the exact BED keep the original, unextended coordinates.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser = clap::value_parser!(u32), help_heading = "Outlier Model")
    )]
    pub blacklist_flank: Option<u32>,

    /// Ignore the gap between paired reads when calculating coverage `[flag]`
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

    /// Optional BED file(s) with positions excluded from coverage estimation and detection `[path]`
    ///
    /// Supply your usual blacklist to detect the sample-specific pileups that remain
    /// outside known problematic regions.
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
            context_weighting: ContextWeighting::Uniform,
            min_context_fragments: DEFAULT_MIN_CONTEXT_FRAGMENTS,
            tail_probability: None,
            tail_probability_multiplier: DEFAULT_TAIL_PROBABILITY_MULTIPLIER,
            target: OutlierTarget::LocalMean,
            model: OutlierModel::Zinb,
            protection_window_size: DEFAULT_PROTECTION_WINDOW_SIZE,
            min_protection_positions: DEFAULT_MIN_PROTECTION_POSITIONS,
            max_outlier_percentage: DEFAULT_MAX_OUTLIER_PERCENTAGE,
            ambiguous_region_action: AmbiguousRegionAction::Preserve,
            plot_preserved_regions: false,
            preserved_plot_context_multiplier: DEFAULT_PRESERVED_PLOT_CONTEXT_MULTIPLIER,
            preserved_plot_min_span: DEFAULT_PRESERVED_PLOT_MIN_SPAN,
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

    /// Set the minimum eligible-base count of each adaptive local fitting context, including zeros.
    pub fn set_bin_size(&mut self, bin_size: u32) {
        self.bin_size = bin_size;
    }

    /// Select how surrounding cores contribute to both local fitting stages.
    pub fn set_context_weighting(&mut self, weighting: ContextWeighting) {
        self.context_weighting = weighting;
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

    /// Enable or disable the optional preserved-region coverage PNGs.
    pub fn set_plot_preserved_regions(&mut self, enabled: bool) {
        self.plot_preserved_regions = enabled;
    }

    /// Set the multiplier used to choose the total preserved-region plot span.
    pub fn set_preserved_plot_context_multiplier(&mut self, multiplier: f64) {
        self.preserved_plot_context_multiplier = multiplier;
    }

    /// Set the minimum reference span of preserved-region plot contexts.
    pub fn set_preserved_plot_min_span(&mut self, span: u32) {
        self.preserved_plot_min_span = span;
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
        self.context_weighting.triangular_fraction()?;
        #[cfg(not(feature = "plotters"))]
        anyhow::ensure!(
            !self.plot_preserved_regions,
            "--plot-preserved-regions requires a build with the plotters feature. Omit this option to write the analysis tables without PNGs"
        );
        anyhow::ensure!(
            self.preserved_plot_context_multiplier.is_finite()
                && self.preserved_plot_context_multiplier >= 1.0,
            "--preserved-plot-context-multiplier must be finite and at least 1"
        );
        anyhow::ensure!(
            self.preserved_plot_min_span > 0,
            "--preserved-plot-min-span must be greater than zero"
        );
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
        push_value(&mut args, "--context-weighting", self.context_weighting);
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
        push_bool(
            &mut args,
            "--plot-preserved-regions",
            self.plot_preserved_regions,
        );
        push_value(
            &mut args,
            "--preserved-plot-context-multiplier",
            self.preserved_plot_context_multiplier,
        );
        push_value(
            &mut args,
            "--preserved-plot-min-span",
            self.preserved_plot_min_span,
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
