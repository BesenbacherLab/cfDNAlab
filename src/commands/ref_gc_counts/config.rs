use crate::commands::cli_common::*;
use std::path::PathBuf;

// TODO: Do we need to add end-offset here, if users use it when calculating gc bias in cfDNA?

/// Count GC fraction per fragment length at a sampled number of starting positions in the reference genome.
/// This 2D count distribution can serve as the expected GC bias in GC correction.
///
/// How: A number (default: 500M) of starting positions are uniformly sampled across the reference
/// genome. For each position, we count the GC fraction for every possible fragment length (default: 30-1000bp).
///
/// Intervals (the possible fragments) with too few ACGT bases after blacklist masking are discarded
/// (so increase `--n-positions` accordingly).
///
/// ## Interpolation
///
/// Some GC fractions are unlikely to see with certain fragment lengths,
/// as only occurence of masked positions will lead to those fractions. Hence, there
/// will be a lot of 0s in the counts. To enable the use of those GC fractions in
/// downstream correction of partly masked fragments, we interpolate the zero-counts (only)
/// using a second-order polynomial with the 3 nearest neighbours on each side.
/// Note: Zeros can still occur at edges.
#[cfg_attr(feature = "cli", derive(clap::Args))]
pub struct RefGCCountsConfig {
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
        clap(long, default_value = "500000000", help_heading = "Core")
    )]
    pub n_positions: usize,

    /// Seed for sampling of start positions `[integer]`
    ///
    /// Use this to reproduce identical reference GC outputs across runs.
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Core"))]
    pub seed: Option<u64>,

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

    /// Number of bases to exclude from each fragment end `[integer]`
    ///
    /// The nucleotides in the cfDNA fragment ends can reflect biological biases (e.g., DNase activity).
    /// This argument allows isolating the GC correction from this signal.
    ///
    /// The default of `10 bp` is based on the GCfix paper by Rahman et al. 2025.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "10",
             value_parser = clap::value_parser!(u8).range(0..20), help_heading="Core"))]
    pub end_offset: u8,

    /// Whether to skip the interpolation of zero-counts `[flag]`
    ///
    /// By default, `0`s are interpolated **independently per fragment length**.
    /// The assumption is that 0s are caused due to the GC content not
    /// being possible to observe with a given fragment length
    /// (e.g., a fragment length of 47 can never achieve a 99% GC).
    /// To avoid errors from this in downstream use, we use polynomial
    /// interpolation based on the neighbourhood of non-zero counts.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub skip_interpolation: bool,

    /// Standard deviation for Gaussian kernel that smoothes raw GC counts for each fragment length `[float]`
    ///
    /// Before converting to discrete GC percentages, we apply smoothing to the raw GC counts separately for each fragment length.
    /// For a fragment length of 150, we thus have counts of fragments with GCs ranging from 0..=150, and smoothing
    /// happens on this scale so the distance between elements are the same for all fragment lengths.
    ///
    /// Note: The same smoothing parameters (sigma and radius) are used for downstream `cfdna gc-bias` calls.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0.55",
             value_parser = clap::value_parser!(f64), help_heading="Smoothing"))]
    pub smoothing_sigma: f64,

    /// Radius of Gaussian kernel that smoothes raw GC counts for each fragment length `[integer]`
    ///
    /// Kernel size is `2 * radius + 1`.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "2",
             value_parser = clap::value_parser!(u8).range(1..10), help_heading="Smoothing"))]
    pub smoothing_radius: u8,

    /// Whether to skip the smoothing of raw GC counts`[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Smoothing"))]
    pub skip_smoothing: bool,
}

impl RefGCCountsConfig {
    pub fn check_smoothing_settings(&self) -> anyhow::Result<()> {
        anyhow::ensure!(
            self.smoothing_sigma > 0.0,
            "--smoothing-sigma must be positive"
        );
        anyhow::ensure!(
            self.smoothing_sigma <= 10.0,
            "--smoothing-sigma must be <= 10.0"
        );
        Ok(())
    }
}
