use std::path::PathBuf;

use crate::commands::cli_common::*;

/// Count fragments per GC fraction and fragment length in a BAM-file.
///
/// Fragment length is defined as `end(reverse) - start(forward)`.
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

    /// 2bit reference file `[path]`
    ///
    /// E.g., "hg38.2bit"
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'r',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub ref_2bit: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub window_assignment: AssignToWindowArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// Masking: Blacklisted positions are set to 'N' in the reference sequence
    /// the GC fraction is calculated from. See the `Minimum ACGT` options
    /// for when to ignore a fragment with too few ACGT (non-'N' and non-blacklisted) bases.
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

    /// Minimum **percentage** of ACGT bases in a fragment after blacklist masking `[integer]`
    ///
    /// Fragments where a lower percentage of bases are ACGT (not blacklisted or 'N') are ignored.
    ///
    /// When both `min_acgt_*` arguments are specified, both thresholds must be met. E.g.,
    /// you may want at least 50% ACGT remaining but also at least 20 bases for a proper
    /// calculation of GC %. For fragments of size 30bp, 50% is only 15bp why the 20bp threshold kicks in.
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
