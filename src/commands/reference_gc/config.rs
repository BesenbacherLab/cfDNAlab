use crate::commands::cli_common::*;
use std::path::PathBuf;

/// Count GC fraction per fragment length at a sampled number of starting positions in the reference genome.
/// This 2D count distribution can serve as the expected GC bias in GC correction.
///
/// How: A number (default: 150M) of starting positions are uniformly sampled across the reference
/// genome. For each position, we count the GC fraction for every possible fragment length (default: 20-1000bp).
///
/// Intervals (the possible fragments) with too few ACGT bases after blacklist masking are discarded
/// (so increase `--n-positions` accordingly).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[cfg_attr(
    feature = "cli",
    clap(
        group = clap::ArgGroup::new("min_acgt")
            .args(&["min_acgt_pct", "min_acgt_count"])
            .multiple(true)))]
pub struct RefGCConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ref_genome: Ref2BitRequiredArgs,

    /// Output directory for results [path]
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'o',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub output_dir: PathBuf,

    /// Number of threads to use (increases RAM usage) [integer]
    ///
    /// Defaults to the minimum of 22 (one thread per chromosome) and
    /// the number of available CPU cores (-1).
    #[cfg_attr(
        feature = "cli",
        clap(short = 't', long, default_value_t = (num_cpus::get()-1).max(1).min(22), help_heading = "Core")
    )]
    pub n_threads: usize,

    /// Number of genomic starting positions to sample [integer]
    ///
    /// The positions are uniformly sampled across the chromosomes
    /// with the GC of each fragment length being counted from
    /// those same starting positions.
    ///
    /// **NOTE**: Sampling is independent of windowing and blacklisting!
    /// The per-length-sum of the output counts may thus be significantly
    /// lower than the specified `n_positions` and different between lengths.
    /// **TIP**: Add 20% extra starting positions than you think you need,
    /// since blacklisting likely removes a big chunk of them.
    #[cfg_attr(
        feature = "cli",
        clap(short = 't', long, default_value = "150000000", help_heading = "Core")
    )]
    pub n_positions: usize,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Optional BED file(s) with blacklisted regions [path]
    ///
    /// Masking: Blacklisted positions are set to 'N' in the reference sequence
    /// the GC fraction is calculated from. See the `Minimum ACGT` options
    /// for when to ignore a kmer with too few ACGT (non-'N' and non-blacklisted) bases.
    #[cfg_attr(
        feature = "cli",
        clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading="Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum **percentage** of ACGT bases in a kmer after blacklist masking [integer]
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

    /// Minimum **count** of ACGT bases in a fragment after blacklist masking [integer]
    ///
    /// Fragments where fewer bases are ACGT (not blacklisted or 'N') are ignored.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20", group = "min_acgt", 
             value_parser = clap::value_parser!(u8).range(0..), help_heading="Minimum ACGT (select 0-2 args)"))]
    pub min_acgt_count: u8,
}
