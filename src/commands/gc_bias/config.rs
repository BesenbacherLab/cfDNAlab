use crate::commands::cli_common::*;
use std::{path::PathBuf, str::FromStr};

#[derive(Default, Clone, Debug)]
pub enum WindowWeightingSchemes {
    Unweighted,
    #[default]
    Coverage,
    Positions,
}

impl WindowWeightingSchemes {
    pub fn as_str(self) -> &'static str {
        match self {
            WindowWeightingSchemes::Unweighted => "unweighted",
            WindowWeightingSchemes::Coverage => "coverage",
            WindowWeightingSchemes::Positions => "positions",
        }
    }
}

impl FromStr for WindowWeightingSchemes {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "unweighted" {
            Ok(WindowWeightingSchemes::Unweighted)
        } else if s == "coverage" {
            Ok(WindowWeightingSchemes::Coverage)
        } else if s == "positions" {
            Ok(WindowWeightingSchemes::Positions)
        } else {
            Err("Use 'unweighted', 'coverage', or 'positions'".into())
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
            .args(&["min_acgt_pct", "min_acgt_count"])
            .multiple(true)))]
#[derive(Clone)]
pub struct GCConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    /// Path to reference GC bias to correct against, calculated with `cfdna reference-gc` `[path]`
    ///
    /// Windowing: When the reference bias is passed in genomic windows, we calculate the
    /// cfDNA GC bias corrections per window and average them (see --average-by).
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, required = true, help_heading = "Core")
    )]
    pub ref_gc: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub window_assignment: AssignToWindowArgs,

    /// How to weight the windows when averaging the correction matrices `[unweighted|coverage|positions]`
    ///
    /// One of:
    ///
    ///  - `"positions"`: Weight windows by how many positions are usable (not blacklisted or `N`).
    ///
    ///  - `"coverage"` (default): Windows are weighted by their average number of observed fragments.
    ///     Compared to a single global window, this approach ensures the local reference bias only
    ///     affects the local correction matrix.
    ///
    ///  - `"unweighted"``: All windows get the same weight in the final correction matrix, no matter how many positions were blacklisted, etc.
    ///  
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "coverage",
            required = true,
            help_heading = "Core"
        )
    )]
    pub window_weighting: WindowWeightingSchemes,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

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

    // TODO: Rethink the effect of this. It will affect smoothing of nearby bins as well.
    // TODO: Should there be a separate global threshold?
    /// Minimum fragment count per combination of GC and fragment length (per window)
    /// before attempting correction `[integer]`
    ///
    /// When a combination of GC content and fragment length has fewer fragments in a given window,
    /// we assign a weight of `1.0`. When averaging, this can reduce extreme corrections from a
    /// few windows where the combination occurs by biasing the correction towards 1.0.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "2",
             value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub min_bin_count: u32,

    /// Minimum fragment count per window to consider it in the correction calculation `[integer]`
    ///
    /// The appropriate threshold depends on the sequencing depth and window size.
    /// On 1 Mb windows, 1000 fragments roughly corresponds to ~0.166x coverage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1000",
             value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub min_window_count: u32,

    /// Minimum percentage of usable positions in a window to consider it in the correction calculation `[integer]`
    ///
    /// Positions are usable if they are not blacklisted or `N`.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20",
             value_parser = clap::value_parser!(u8).range(0..101), help_heading="Filtering"))]
    pub min_usable_positions_pct: u8,

    /// Minimum **percentage** of ACGT bases in a fragment after blacklist masking `[integer]`
    ///
    /// Fragments where a lower percentage of bases are ACGT (not blacklisted or 'N') are ignored.
    ///
    /// When both `min_acgt_*` arguments are specified, both thresholds must be met. E.g.,
    /// you may want at least 50% ACGT remaining but also at least 20 bases for a proper
    /// calculation of GC %. For fragments of size 30bp, 50% is only 15bp, so the 20bp
    /// absolute threshold kicks in.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "90", group = "min_acgt", 
             value_parser = clap::value_parser!(u8).range(0..101), help_heading="Minimum ACGT (select 0-2 args)"))]
    pub min_acgt_pct: u8,

    /// Minimum **count** of ACGT bases in a fragment after blacklist masking `[integer]`
    ///
    /// Fragments where fewer bases are ACGT (not blacklisted or 'N') are ignored.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20", group = "min_acgt", 
             value_parser = clap::value_parser!(u8).range(0..), help_heading="Minimum ACGT (select 0-2 args)"))]
    pub min_acgt_count: u8,
}
