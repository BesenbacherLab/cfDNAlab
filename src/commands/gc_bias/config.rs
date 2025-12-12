use crate::commands::cli_common::*;
use crate::commands::gc_bias::outliers::{OutlierAction, OutlierRule, OutlierScope};
use anyhow::{Result, anyhow};
use std::{path::PathBuf, str::FromStr};

#[derive(Default, Clone, Debug)]
pub enum WindowWeightingSchemes {
    Equal,
    Coverage,
    #[default]
    ValidPositions,
}

impl WindowWeightingSchemes {
    pub fn as_str(self) -> &'static str {
        match self {
            WindowWeightingSchemes::Equal => "equal",
            WindowWeightingSchemes::Coverage => "coverage",
            WindowWeightingSchemes::ValidPositions => "valid-positions",
        }
    }
}

impl FromStr for WindowWeightingSchemes {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "equal" {
            Ok(WindowWeightingSchemes::Equal)
        } else if s == "coverage" {
            Ok(WindowWeightingSchemes::Coverage)
        } else if s == "valid-positions" {
            Ok(WindowWeightingSchemes::ValidPositions)
        } else {
            Err("Use 'equal', 'coverage', or 'valid-positions'".into())
        }
    }
}

#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutlierMethodArg {
    None,
    Quantile,
    Iqr,
    Stddev,
    Mad,
}

impl Default for OutlierMethodArg {
    fn default() -> Self {
        OutlierMethodArg::Iqr
    }
}

#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutlierScopeArg {
    PerLength,
    Global,
}

impl Default for OutlierScopeArg {
    fn default() -> Self {
        OutlierScopeArg::PerLength
    }
}

// TODO: Try excluding the first N bases (both ends) from GC fraction calculation to avoid correcting "biochemical cut bias" - the bias we care about is "regional bias"
// Perhaps do an "end-proximal base composition (p≈1–10) bias" experiment to show how many bases to cut off in ends

/// Calculate a multiplicative GC correction matrix based on the GC fraction and length of fragments in a BAM-file.
///
/// The observed distribution of cfDNA fragments is corrected to a precomputed reference bias.
///
/// Requirements: Please precompute the reference GC bias with `cfdna reference-gc`.
/// This file can be reused for all samples (aligned to the same assembly).
///
/// The most extreme GC bins get corrections of `1.0` to avoid extreme corrections due to sparsity.
///
/// Combinations of GC fractions and fragment lengths that are either theoretically unobservable
/// or *very* rarely observed in the *reference genome* are interpolated from surrounding counts.
/// Other combinations with zero counts in the *cfDNA* remains zero in the correction matrix.
/// The final correction matrix thus works for all possible GC x Length combinations.
///
/// Fragment length is defined as `end(reverse) - start(forward)`.
///
/// ## Windowing
///
/// Technical GC bias is assumed to be a "global" bias. To control how each region of the genome
/// (which may have amplified or reduced coverage) contributes to the calculation of this global bias,
/// we can calculate the bias in genomic windows and combine them via (weighted) averaging.
///
/// The windows are taken from the reference GC bias file from `cfdna reference-gc`.
///
/// The output of `cfdna gc-bias` is always a 2D correction matrix.
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
pub struct GCConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ref_genome: Ref2BitRequiredArgs,

    /// Path to directory with reference GC bias files to correct against `[path]`
    ///
    /// Precompute with `cfdna reference-gc`. The directory must include all files
    /// created by `cfdna reference-gc`.
    ///
    /// Windowing: When the reference bias is passed in genomic windows, we calculate the
    /// cfDNA GC bias corrections per window and average them (see `--window-weighting`).
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, required = true, help_heading = "Core")
    )]
    pub ref_gc_dir: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub window_assignment: AssignToWindowArgs,

    /// How to weight the windows when averaging them `[coverage|valid-positions|equal]`
    ///
    /// One of:
    ///
    ///  - `"coverage"`: Windows are weighted by their average number of observed fragments.
    ///    Compared to a single global window, this approach weights the local reference bias the
    ///    same as the local cfDNA bias in the global biases. *Technically*, only the reference count
    ///    distribution is reweighted, as the cfDNA counts already reflect the coverage.
    ///
    ///  - `"valid-positions"`: Weight windows by how many positions are usable (not blacklisted or `N`).
    ///
    ///  - `"equal"`: All windows get the same weight in the final correction matrix,
    ///    no matter how many positions were blacklisted, etc.
    ///
    /// **NOTE**: Only specify this argument when windows exist.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_enum, default_value = "valid-positions", help_heading = "Core")
    )]
    pub window_weighting: WindowWeightingSchemes,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    /// Minimum percentage of counts to have in each length bin `[float]`
    ///
    /// Greater than 0, lower than 100. Default is 0.5% (i.e., a max. of 200 bins).
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0.5", value_parser = parse_percentage_within_0_100_f32, help_heading="Binning"))]
    pub min_length_bin_mass: f32,

    /// Minimum number of fragment lengths per fragment length bin `[float]`
    ///
    /// Reduces sparsity-related issues in ultra low-coverage samples.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "3", value_parser = clap::value_parser!(u8).range(1..100), help_heading="Binning"))]
    pub min_length_bin_width: u8,

    /// Minimum percentage of counts to have in each GC contents bin `[float]`
    ///
    /// Greater than 0, lower than 100. Default is 1% (i.e., a max. of 100 bins).
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1.0", value_parser = parse_percentage_within_0_100_f32, help_heading="Binning"))]
    pub min_gc_bin_mass: f32,

    /// Number of extreme GC bins (`--min_gc_bin_mass`) from each side to interpolate from neighbouring corrections `[float]`
    ///
    /// The most extreme GC fractions are very sparsely observed. This can lead to extreme corrections.
    /// Set the number of bins from each side where we interpolate a correction based on the neighbouring corrections.
    /// The default of 1 should be fine but this can be tuned via visualization of the created
    /// correction matrix and intermediate files (`--save-intermediates`).
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1", value_parser = clap::value_parser!(u8).range(0..10), help_heading="Binning"))]
    pub num_extreme_gc_bins: u8,

    /// Number of the **shortest** fragment length bins (`--min_length_bin_mass`) to interpolate from neighbouring corrections `[float]`
    ///
    /// The shortest fragment lengths can be very sparsely observed. This can lead to extreme corrections.
    /// Set the number of short-length bins where we interpolate a correction based on the neighbouring corrections.
    /// With the default minimum fragment length setting in `cfdna reference-gc` (30bp),
    /// the default of 1 should be fine. This can be tuned via visualization of the created
    /// correction matrix and intermediate files (`--save-intermediates`).
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1", value_parser = clap::value_parser!(u8).range(0..10), help_heading="Binning"))]
    pub num_short_length_bins: u8,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// Masking: Blacklisted positions are set to 'N' in the reference sequence
    /// the GC fraction is calculated from. See the `Minimum ACGT` options
    /// for when to ignore a fragment with too few ACGT (non-'N' and non-blacklisted) bases.
    ///
    /// NOTE: Ensure the same positions were blacklisted when calculating the reference bias (`cfdna reference-gc`).
    #[cfg_attr(
        feature = "cli",
        clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading="Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

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

    /// Minimum percentage of ACGT positions in a **window** to consider it in the bias estimation `[integer]`
    ///
    /// If you believe windows that are mostly blacklisted may be too noisy in their
    /// remaining positions, use this to threshold to remove them from the analysis.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "10",
             value_parser = clap::value_parser!(u8).range(0..101), help_heading="Minimum ACGT"))]
    pub min_window_acgt_pct: u8,

    /// Handle extreme correction factors to avoid unstable weights `[string]`
    ///
    /// Options:
    ///
    /// - `none`: Disable outlier handling.
    ///
    /// - `quantile`: Clamp using `--outlier-quantiles` (one symmetric value or two explicit values).
    ///
    /// - `iqr`, `stddev`, `mad`: Use the corresponding rule with multiplier `--outlier-k`.
    ///
    /// **NOTE**: After outlier detection, correction values are further clipped at `[0.1, 10.0]`.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "iqr", value_enum, help_heading = "Outliers")
    )]
    pub outlier_method: OutlierMethodArg,

    /// Whether to detect outliers per fragment length or across the full matrix `[string]`
    ///
    /// - `per-length`: Detect separately per fragment length.
    ///
    /// - `global`: Detect from the full correction matrix.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "global", value_enum, help_heading = "Outliers")
    )]
    pub outlier_scope: OutlierScopeArg,

    /// Quantiles for `quantile` outlier detection `[float or float,float]`
    ///
    /// Used when `--outlier-method quantile`. Provide one value to apply symmetrically (`q` -> lower=`q`, upper=`1-q`)
    /// or two values for explicit `lower,upper`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = parse_quantile_0_1,
            num_args = 1..=2,
            default_values_t = [0.03_f32, 0.97_f32],
            help_heading = "Outliers"
        )
    )]
    pub outlier_quantiles: Vec<f32>,

    /// Multiplier `k` for `iqr`, `stddev`, or `mad` outlier detection `[float]`
    ///
    /// Used when `--outlier-method` is one of `iqr`, `stddev`, or `mad`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "3",
            value_parser = clap::value_parser!(f32),
            help_heading = "Outliers"
        )
    )]
    pub outlier_k: f32,

    /// Whether to save key intermediate files for inspecting the correction process `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub save_intermediates: bool,
}

impl GCConfig {
    pub fn new(
        ioc: IOCArgs,
        ref_2bit: PathBuf,
        ref_gc_dir: PathBuf,
        chromosomes: ChromosomeArgs,
    ) -> Self {
        Self {
            ioc,
            ref_genome: Ref2BitRequiredArgs { ref_2bit },
            ref_gc_dir,
            window_assignment: AssignToWindowArgs::default(),
            window_weighting: WindowWeightingSchemes::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            blacklist: None,
            min_mapq: 30,
            require_proper_pair: false,
            min_gc_bin_mass: 1.0,
            min_length_bin_mass: 1.0,
            min_length_bin_width: 3,
            num_extreme_gc_bins: 1,
            num_short_length_bins: 1,
            min_window_acgt_pct: 10,
            outlier_method: OutlierMethodArg::Iqr,
            outlier_scope: OutlierScopeArg::PerLength,
            outlier_quantiles: vec![0.03, 0.97],
            outlier_k: 8.0,
            save_intermediates: false,
        }
    }

    pub fn set_ioc(&mut self, ioc: IOCArgs) {
        self.ioc = ioc;
    }

    pub fn set_ref_gc_dir(&mut self, ref_gc_dir: PathBuf) {
        self.ref_gc_dir = ref_gc_dir;
    }

    pub fn set_window_assignment(&mut self, assignment: AssignToWindowArgs) {
        self.window_assignment = assignment;
    }

    pub fn set_window_weighting(&mut self, weighting: WindowWeightingSchemes) {
        self.window_weighting = weighting;
    }

    pub fn set_chromosomes(&mut self, chromosomes: ChromosomeArgs) {
        self.chromosomes = chromosomes;
    }

    pub fn set_scaling_factors(&mut self, scaling_factors: Option<PathBuf>) {
        self.scale_genome.scaling_factors = scaling_factors;
    }

    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_min_length_bin_mass(&mut self, min_length_bin_mass: f32) {
        self.min_length_bin_mass = min_length_bin_mass;
    }

    pub fn set_min_length_bin_width(&mut self, min_length_bin_width: u8) {
        self.min_length_bin_width = min_length_bin_width;
    }

    pub fn set_min_gc_bin_mass(&mut self, min_gc_bin_mass: f32) {
        self.min_gc_bin_mass = min_gc_bin_mass;
    }

    pub fn set_num_extreme_gc_bins(&mut self, num_extreme_gc_bins: u8) {
        self.num_extreme_gc_bins = num_extreme_gc_bins;
    }

    pub fn set_num_short_length_bins(&mut self, num_short_length_bins: u8) {
        self.num_short_length_bins = num_short_length_bins;
    }

    pub fn set_min_window_acgt_pct(&mut self, pct: u8) {
        self.min_window_acgt_pct = pct;
    }

    pub fn set_save_intermediates(&mut self, save_intermediates: bool) {
        self.save_intermediates = save_intermediates;
    }

    pub fn outlier_settings(&self) -> Result<(OutlierRule, OutlierAction, OutlierScope)> {
        let rule = match self.outlier_method {
            OutlierMethodArg::None => OutlierRule::None,
            OutlierMethodArg::Quantile => {
                let (lower, upper) = resolve_outlier_quantiles(&self.outlier_quantiles)?;
                OutlierRule::Quantile { lower, upper }
            }
            OutlierMethodArg::Iqr => {
                if self.outlier_k <= 0.0 {
                    return Err(anyhow!("outlier-k must be > 0 for IQR"));
                }
                OutlierRule::TukeyIqr { k: self.outlier_k }
            }
            OutlierMethodArg::Stddev => {
                if self.outlier_k <= 0.0 {
                    return Err(anyhow!("outlier-k must be > 0 for stddev"));
                }
                OutlierRule::StdDev { k: self.outlier_k }
            }
            OutlierMethodArg::Mad => {
                if self.outlier_k <= 0.0 {
                    return Err(anyhow!("outlier-k must be > 0 for mad"));
                }
                OutlierRule::Mad { k: self.outlier_k }
            }
        };

        let action = OutlierAction::Winsorize;

        let scope = match self.outlier_scope {
            OutlierScopeArg::PerLength => OutlierScope::PerLength,
            OutlierScopeArg::Global => OutlierScope::Global,
        };

        Ok((rule, action, scope))
    }
}

#[cfg_attr(not(feature = "cli"), allow(dead_code))]
fn parse_percentage_within_0_100_f32(input: &str) -> Result<f32, String> {
    let value: f32 = input
        .parse()
        .map_err(|err: std::num::ParseFloatError| err.to_string())?;
    if value <= 0.0 {
        Err("value must be > 0".into())
    } else if value >= 100.0 {
        Err("value must be < 100".into())
    } else {
        Ok(value)
    }
}

fn parse_quantile_0_1(input: &str) -> Result<f32, String> {
    let value: f32 = input
        .parse()
        .map_err(|err: std::num::ParseFloatError| err.to_string())?;
    if !(0.0..=1.0).contains(&value) {
        Err("value must be between 0 and 1".into())
    } else {
        Ok(value)
    }
}

fn resolve_outlier_quantiles(vals: &[f32]) -> Result<(f32, f32)> {
    if vals.is_empty() {
        return Err(anyhow!(
            "provide at least one quantile for --outlier-method quantile"
        ));
    }
    if vals.len() == 1 {
        let q = vals[0];
        if q >= 0.5 {
            return Err(anyhow!(
                "single quantile must be < 0.5 to allow symmetric bounds (got {})",
                q
            ));
        }
        return Ok((q, 1.0 - q));
    }
    let lower = vals[0];
    let upper = vals[1];
    if lower >= upper {
        return Err(anyhow!(
            "outlier lower quantile must be < upper quantile (got {} >= {})",
            lower,
            upper
        ));
    }
    Ok((lower, upper))
}
