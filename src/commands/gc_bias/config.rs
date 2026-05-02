use crate::commands::cli_common::*;
use crate::commands::gc_bias::outliers::{OutlierAction, OutlierRule, OutlierScope};
use anyhow::{Result, anyhow};
use std::path::PathBuf;

// Central defaults to keep CLI and programmatic creation in sync
pub const DEFAULT_TILE_SIZE: u32 = 10_000_000;
pub const DEFAULT_MIN_LENGTH_BIN_MASS: f32 = 0.5;
pub const DEFAULT_MIN_LENGTH_BIN_WIDTH: u8 = 3;
pub const DEFAULT_MIN_GC_BIN_MASS: f32 = 1.0;
pub const DEFAULT_NUM_EXTREME_GC_BINS: u8 = 1;
pub const DEFAULT_NUM_SHORT_LENGTH_BINS: u8 = 1;
pub const DEFAULT_MIN_WINDOW_ACGT_PCT: u8 = 10;
pub const DEFAULT_MIN_MAPQ: u8 = 30;
pub const DEFAULT_OUTLIER_K: f32 = 3.0;
pub const DEFAULT_OUTLIER_QUANTILES: [f32; 2] = [0.03, 0.97];
pub const DEFAULT_OUTLIER_METHOD: OutlierMethodArg = OutlierMethodArg::Iqr;
pub const DEFAULT_OUTLIER_SCOPE: OutlierScopeArg = OutlierScopeArg::Global;

#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutlierMethodArg {
    None,
    Quantile,
    #[default]
    Iqr,
    Stddev,
    Mad,
}

#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutlierScopeArg {
    PerLength,
    #[default]
    Global,
}

/// Calculate a multiplicative GC correction matrix based on the GC fraction and length of fragments in a BAM-file.
///
/// The observed distribution of cfDNA fragments is corrected to a precomputed reference bias.
///
/// Requirements: Please precompute the reference GC bias with `cfdna ref-gc-bias`.
/// This file can be reused for all samples aligned to the same assembly.
///
/// **NOTE**: This command is highly flexible, enabling experimentation. The default values have been
/// tuned and should be useful in most use cases. Start with the example below.
///
/// ## Interpolations
///
/// The most extreme GC and shortest-length bins get interpolated corrections based on neighbours
/// to avoid extreme corrections due to sparsity.
///
/// The combinations of GC fractions and fragment lengths that are either theoretically unobservable
/// or *very* rarely observed in the **reference genome** are interpolated from surrounding counts.
/// Other combinations with post-smoothing zero counts in the *cfDNA* remains zero in the correction matrix.
/// The final correction matrix thus works for all possible GC x Length combinations.
///
/// ## Fragment length definition
///
/// **Paired-end**: `end(reverse) - start(forward)`.
///
/// **Unpaired** where each read is a fragment: `end(read) - start(read)`.
///
/// The utilized fragment length range is inherited from the
/// `--ref-gc-file` to ensure consistency.
///
/// ## Windowing
///
/// Technical GC bias is assumed to be a "global" bias. To control how each region of the genome
/// (which may have amplified/reduced coverage) contributes to the calculation of this global bias,
/// we can calculate the bias in genomic windows and combine them via weighted averaging:
/// The counts of each window are divided by their window-mean and scaled by the number of
/// valid ACGT positions in the window. The windows are then averaged.
///
/// ## Example
///
/// ```bash
///
/// cfdna gc-bias --bam {BAM_FILE} --output-dir {PATH}/gc_bias \
///
///   --ref-2bit {PATH}/hg38.2bit \ # Or some other assembly
///
///   --ref-gc-file {REFERENCE_GC_FILE} \
///
///   --blacklist {PATH}/encode_blacklist.bed # Or some other blacklist(s)
///
/// ```
///
/// Besides these arguments, the default values should work in most cases.
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
pub struct GCConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.gc_bias_correction.npz`
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, help_heading = "Core")
    )]
    pub output_prefix: String,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ref_genome: Ref2BitRequiredArgs,

    /// Path to file with reference GC bias to correct against `[path]`
    ///
    /// Precompute with `cfdna ref-gc-bias`. The file is either named
    /// `ref_gc_package.npz` or `<prefix>.ref_gc_package.npz`.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, required = true, help_heading = "Core")
    )]
    pub ref_gc_file: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: GCWindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub window_assignment: AssignToWindowArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_TILE_SIZE,
            value_parser = clap::value_parser!(u32).range(1000000..),
            help_heading = "Core"
        )
    )]
    pub tile_size: u32,

    /// Minimum percentage of counts to have in each length bin `[float]`
    ///
    /// Greater than 0, lower than 100. Default is 0.5% (i.e., a max. of 200 bins).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_MIN_LENGTH_BIN_MASS,
            value_parser = parse_percentage_within_0_100_f32,
            help_heading = "Binning"
        )
    )]
    pub min_length_bin_mass: f32,

    /// Minimum number of fragment lengths per fragment length bin `[float]`
    ///
    /// Reduces sparsity-related issues in ultra low-coverage samples.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_MIN_LENGTH_BIN_WIDTH,
            value_parser = clap::value_parser!(u8).range(1..100),
            help_heading = "Binning"
        )
    )]
    pub min_length_bin_width: u8,

    /// Minimum percentage of counts to have in each GC contents bin `[float]`
    ///
    /// Greater than 0, lower than 100. Default is 1% (i.e., a max. of 100 bins).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_MIN_GC_BIN_MASS,
            value_parser = parse_percentage_within_0_100_f32,
            help_heading = "Binning"
        )
    )]
    pub min_gc_bin_mass: f32,

    /// Number of extreme GC bins (`--min_gc_bin_mass`) from each side to interpolate from neighbouring corrections `[integer]`
    ///
    /// The most extreme GC fractions are very sparsely observed. This can lead to extreme corrections.
    /// Set the number of bins from each side where we interpolate a correction based on the neighbouring corrections.
    /// The default of 1 should be fine but this can be tuned via visualization of the created
    /// correction matrix and intermediate files (`--save-intermediates`).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_NUM_EXTREME_GC_BINS,
            value_parser = clap::value_parser!(u8).range(0..10),
            help_heading = "Binning"
        )
    )]
    pub num_extreme_gc_bins: u8,

    /// Number of the **shortest** fragment length bins (`--min_length_bin_mass`) to interpolate from neighbouring corrections `[integer]`
    ///
    /// The shortest fragment lengths can be very sparsely observed. This can lead to extreme corrections.
    /// Set the number of short-length bins where we interpolate a correction based on the neighbouring corrections.
    /// With the default minimum fragment length setting in `cfdna ref-gc-bias` (30bp),
    /// the default of 1 should be fine. This can be tuned via visualization of the created
    /// correction matrix and intermediate files (`--save-intermediates`).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_NUM_SHORT_LENGTH_BINS,
            value_parser = clap::value_parser!(u8).range(0..10),
            help_heading = "Binning"
        )
    )]
    pub num_short_length_bins: u8,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// Masking: Blacklisted positions are set to 'N' in the reference sequence
    /// the GC fraction is calculated from. See the `Minimum ACGT` options
    /// for when to ignore a fragment with too few ACGT (non-'N' and non-blacklisted) bases.
    ///
    /// NOTE: Ensure the same positions were blacklisted when calculating the reference bias (`cfdna ref-gc-bias`).
    #[cfg_attr(
        feature = "cli",
        clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading="Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "mq",
            default_value_t = DEFAULT_MIN_MAPQ,
            value_parser = clap::value_parser!(u8).range(0..),
            help_heading = "Filtering"
        )
    )]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    ///
    /// This is **NOT** recommended by default as it trims the tails of the length distribution.
    ///
    /// Note, that we only keep inward-directed fragments within the specified length range, so
    /// there's no real need for proper-pair filtering.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Minimum percentage of ACGT positions in a **window** to consider it in the bias estimation `[integer]`
    ///
    /// If you believe windows that are mostly blacklisted may be too noisy in their
    /// remaining positions, use this to threshold to remove them from the analysis.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_MIN_WINDOW_ACGT_PCT,
            value_parser = clap::value_parser!(u8).range(0..101),
            help_heading = "Minimum ACGT"
        )
    )]
    pub min_window_acgt_pct: u8,

    /// Handle extreme GC-bias values to avoid unstable weights `[string]`
    ///
    /// Options:
    ///
    /// - `none`: Disable outlier handling.
    ///
    /// - `quantile`: Clamp using `--outlier-quantiles` (one symmetric value or two explicit values).
    ///
    /// - `iqr`, `stddev`, `mad`: Use the corresponding rule with multiplier `--outlier-k`.
    ///
    /// **NOTE**: After outlier detection, extreme GC-bias values are clipped at `[0.1, 10.0]`
    /// before the final scaling steps.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_OUTLIER_METHOD, value_enum, help_heading = "Outliers")
    )]
    pub outlier_method: OutlierMethodArg,

    /// Whether to detect outliers per fragment length or across the full matrix `[string]`
    ///
    /// - `per-length`: Detect separately per fragment length.
    ///
    /// - `global`: Detect from the full correction matrix.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_OUTLIER_SCOPE, value_enum, help_heading = "Outliers")
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
            default_values_t = DEFAULT_OUTLIER_QUANTILES,
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
            default_value_t = DEFAULT_OUTLIER_K,
            value_parser = clap::value_parser!(f32),
            help_heading = "Outliers"
        )
    )]
    pub outlier_k: f32,

    /// Whether to save key intermediate files for inspecting the correction process `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub save_intermediates: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl GCConfig {
    pub fn new(
        ioc: IOCArgs,
        ref_2bit: PathBuf,
        ref_gc_file: PathBuf,
        chromosomes: ChromosomeArgs,
    ) -> Self {
        Self {
            ioc,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix: String::new(),
            ref_genome: Ref2BitRequiredArgs { ref_2bit },
            ref_gc_file,
            windows: GCWindowsArgs::default(),
            window_assignment: AssignToWindowArgs::default(),
            chromosomes,
            tile_size: DEFAULT_TILE_SIZE,
            blacklist: None,
            min_mapq: DEFAULT_MIN_MAPQ,
            require_proper_pair: false,
            min_gc_bin_mass: DEFAULT_MIN_GC_BIN_MASS,
            min_length_bin_mass: DEFAULT_MIN_LENGTH_BIN_MASS,
            min_length_bin_width: DEFAULT_MIN_LENGTH_BIN_WIDTH,
            num_extreme_gc_bins: DEFAULT_NUM_EXTREME_GC_BINS,
            num_short_length_bins: DEFAULT_NUM_SHORT_LENGTH_BINS,
            min_window_acgt_pct: DEFAULT_MIN_WINDOW_ACGT_PCT,
            outlier_method: DEFAULT_OUTLIER_METHOD,
            outlier_scope: DEFAULT_OUTLIER_SCOPE,
            outlier_quantiles: DEFAULT_OUTLIER_QUANTILES.to_vec(),
            outlier_k: DEFAULT_OUTLIER_K,
            save_intermediates: false,
            logging: LoggingArgs::default(),
        }
    }

    pub fn set_ioc(&mut self, ioc: IOCArgs) {
        self.ioc = ioc;
    }

    pub fn set_output_prefix(&mut self, output_prefix: String) {
        self.output_prefix = output_prefix;
    }

    pub fn set_ref_gc_file(&mut self, ref_gc_file: PathBuf) {
        self.ref_gc_file = ref_gc_file;
    }

    pub fn set_windows(&mut self, windows: GCWindowsArgs) {
        self.windows = windows;
    }

    pub fn set_window_assignment(&mut self, assignment: AssignToWindowArgs) {
        self.window_assignment = assignment;
    }

    pub fn set_chromosomes(&mut self, chromosomes: ChromosomeArgs) {
        self.chromosomes = chromosomes;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
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

#[cfg_attr(not(feature = "cli"), allow(dead_code))]
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
