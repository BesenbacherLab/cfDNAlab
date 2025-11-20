use crate::commands::cli_common::*;
use std::{path::PathBuf, str::FromStr};

#[derive(Default, Clone, Debug)]
pub enum WindowWeightingSchemes {
    Equal,
    #[default]
    Coverage,
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

// TODO: Try excluding the first N bases (both ends) from GC fraction calculation to avoid correcting "biochemical cut bias" - the bias we care about is "regional bias"
// Perhaps do an "end-proximal base composition (p≈1–10) bias" experiment to show how many bases to cut off in ends

/// Count fragments per GC fraction and fragment length in a BAM-file.
///
/// Fragment length is defined as `end(reverse) - start(forward)`.
///
/// Requirements: Please precompute the reference GC bias with `cfdna reference-gc`.
/// This file can be reused for all samples (aligned to the same assembly).
///
/// ## Windowing
///
/// Technical GC bias is assumed to be a "global" bias. To control how each region of the genome
/// (which may have amplified or reduced coverage) contributes to the calculation of this global bias,
/// we can calculate the bias in genomic windows and combine them via (weighted) averaging.
///
/// The windows are taken from the reference GC bias file from `cfdna reference-gc`.
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
#[cfg_attr(
    feature = "cli",
    clap(
        group = clap::ArgGroup::new("min_acgt")
            .args(&["min_fragment_acgt_pct", "min_fragment_acgt_count"])
            .multiple(true)))]
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

    /// How to weight the windows when averaging the windows `[valid-positions|coverage|equal]`
    ///
    /// One of:
    ///
    ///  - `"coverage"` (default): Windows are weighted by their average number of observed fragments.
    ///    Compared to a single global window, this approach weights the local reference bias the
    ///    same as the local cfDNA bias in the global biases. *Technically*, only the reference count
    ///    distribution is reweighted, as the cfDNA counts already reflect the coverage.
    ///
    ///  - `"valid-positions"`: Weight windows by how many positions are usable (not blacklisted or `N`).
    ///    This gives equal weight to all valid positions in the genome covered by the windows
    ///    (assuming no overlap between windows).
    ///
    ///  - `"equal"`: All windows get the same weight in the final correction matrix,
    ///    no matter how many positions were blacklisted, etc.
    ///
    /// **NOTE**: Only specify this argument when windows exist.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_enum, default_value = "coverage", help_heading = "Core")
    )]
    pub window_weighting: WindowWeightingSchemes,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    /// Minimum percentage of counts to have in each length bin `[float]`
    ///
    /// Greater than 0, lower than 100.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1.0", value_parser = parse_percentage_within_0_100_f32, help_heading="Binning"))]
    pub min_length_bin_mass: f32,

    /// Minimum percentage of counts to have in each GC contents bin `[float]`
    ///
    /// Greater than 0, lower than 100.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1.0", value_parser = parse_percentage_within_0_100_f32, help_heading="Binning"))]
    pub min_gc_bin_mass: f32,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// Masking: Blacklisted positions are set to 'N' in the reference sequence
    /// the GC fraction is calculated from. See the `Minimum ACGT` options
    /// for when to ignore a fragment with too few ACGT (non-'N' and non-blacklisted) bases.
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

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum GC % to consider `[integer]`
    ///
    /// Fragments with lower GC % are ignored.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0", 
             value_parser = clap::value_parser!(u8).range(0..100), help_heading="Filtering"))]
    pub gc_min_pct: u8,

    /// Maximum GC % to consider `[integer]`
    ///
    /// Fragments with higher GC % are ignored.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "100", 
             value_parser = clap::value_parser!(u8).range(0..101), help_heading="Filtering"))]
    pub gc_max_pct: u8,

    // TODO: Base this on the original GC paper. Look it up. Perhaps add a reference.
    /// Number of bases to exclude from each fragment end `[integer]`
    ///
    /// The nucleotides in the fragment ends can reflect biological biases (e.g., DNase activity).
    /// This argument allows isolating the GC correction from this signal.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0",
             value_parser = clap::value_parser!(u8).range(0..20), help_heading="Filtering"))]
    pub end_offset: u8,

    /// Minimum percentage of ACGT positions in a **window** to consider it in the bias estimation `[integer]`
    ///
    /// If you believe windows that are mostly blacklisted may be too noisy in their
    /// remaining positions, use this to threshold to remove them from the analysis.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "10",
             value_parser = clap::value_parser!(u8).range(0..101), help_heading="Filtering"))]
    pub min_window_acgt_pct: u8,

    /// Minimum **percentage** of ACGT bases in a fragment after blacklist masking and end offsets `[integer]`
    ///
    /// Fragments where a lower percentage of bases are ACGT (not blacklisted or 'N') are ignored.
    /// When specifying `--end-offset`, the ends are excluded before this calculation.
    ///
    /// When both `min_acgt_*` arguments are specified, both thresholds must be met. E.g.,
    /// you may want at least 50% ACGT remaining but also at least 20 bases for a proper
    /// calculation of GC %. For fragments of size 30bp, 50% is only 15bp, so the 20bp
    /// absolute threshold kicks in.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "90", group = "min_acgt", 
             value_parser = clap::value_parser!(u8).range(0..101), help_heading="Minimum ACGT (select 0-2 args)"))]
    pub min_fragment_acgt_pct: u8,

    /// Minimum **count** of ACGT bases in a fragment after blacklist masking and end offsets `[integer]`
    ///
    /// Fragments where fewer bases are ACGT (not blacklisted or 'N') are ignored.
    /// When specifying `--end-offset`, the ends are excluded before this calculation.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20", group = "min_acgt", 
             value_parser = clap::value_parser!(u8).range(0..), help_heading="Minimum ACGT (select 0-2 args)"))]
    pub min_fragment_acgt_count: u8,
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
            fragment_lengths: FragmentLengthArgs {
                min_fragment_length: 20,
                max_fragment_length: 1000,
            },
            gc_min_pct: 0,
            gc_max_pct: 100,
            end_offset: 0,
            min_gc_bin_mass: 1.0,
            min_length_bin_mass: 1.0,
            min_window_acgt_pct: 10,
            min_fragment_acgt_pct: 90,
            min_fragment_acgt_count: 20,
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

    pub fn set_fragment_lengths(&mut self, fragment_lengths: FragmentLengthArgs) {
        self.fragment_lengths = fragment_lengths;
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_gc_min_pct(&mut self, gc_min_pct: u8) {
        self.gc_min_pct = gc_min_pct;
    }

    pub fn set_gc_max_pct(&mut self, gc_max_pct: u8) {
        self.gc_max_pct = gc_max_pct;
    }

    pub fn set_end_offset(&mut self, end_offset: u8) {
        self.end_offset = end_offset;
    }

    pub fn set_min_length_bin_mass(&mut self, min_length_bin_mass: f32) {
        self.min_length_bin_mass = min_length_bin_mass;
    }

    pub fn set_min_gc_bin_mass(&mut self, min_gc_bin_mass: f32) {
        self.min_gc_bin_mass = min_gc_bin_mass;
    }

    pub fn set_min_window_acgt_pct(&mut self, pct: u8) {
        self.min_window_acgt_pct = pct;
    }

    pub fn set_min_fragment_acgt_pct(&mut self, pct: u8) {
        self.min_fragment_acgt_pct = pct;
    }

    pub fn set_min_fragment_acgt_count(&mut self, count: u8) {
        self.min_fragment_acgt_count = count;
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
