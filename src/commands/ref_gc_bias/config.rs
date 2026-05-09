use crate::commands::cli_common::*;
use std::path::PathBuf;

/// Build a reference GC bias table for cfDNA correction.
///
/// Samples approximately `n_positions` across all chromosomes and counts GC for every fragment length in range
/// (optionally trimmed in ends). Creates one genome-wide GC-by-length table that
/// downstream GC bias correction uses as the expected bias. If you provide a BED file via `--by-bed`,
/// overlapping intervals are merged and counting is limited to those bases. Problematic regions
/// can be excluded via a blacklist. Otherwise, the full genome is used.
///
/// After counting, the table is smoothed length-wise and converted to GC percentages.
/// A support mask flags bins with too few counts per megabase (including theoretically unobservable
/// GC-by-length combinations), and the sparse bins are interpolated using neighbours.
#[cfg_attr(feature = "cli", derive(clap::Args))]
pub struct RefGCBiasConfig {
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

    /// Optional prefix for the output file (e.g., a reference genome name) `[string]`
    ///
    /// Leave empty to write the filename without a leading prefix.
    ///
    /// E.g., to allow storing packages for multiple reference genomes in the same directory.
    ///
    /// Produces the file as:
    ///   `<prefix>.ref_gc_package.npz`
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Number of threads to use (increases RAM usage) [integer]
    ///
    /// Defaults to the number of available CPU cores (-1).
    #[cfg_attr(
        feature = "cli",
        clap(short = 't', long, default_value_t = (num_cpus::get()-1).max(1), help_heading = "Core")
    )]
    pub n_threads: usize,

    /// Number of genomic starting positions to sample [integer]
    ///
    /// The positions are uniformly sampled across the chromosomes
    /// with the GC of each fragment length being counted from
    /// those same starting positions.
    ///
    /// **NOTE**: `--n-positions` is an approximate sampling target, not an exact quota.
    /// Sampling is independent of windowing and blacklisting and the per-length-sum
    /// of the output counts may thus be significantly lower than the specified
    /// `n_positions` and different between lengths.
    ///
    /// **TIP**: Add 20% extra starting positions than you think you need,
    /// since blacklisting likely removes a big chunk of them.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "500000000",
            value_parser = parse_positive_usize,
            help_heading = "Core"
        )
    )]
    pub n_positions: usize,

    /// Seed for sampling of start positions `[integer]`
    ///
    /// Use this to reproduce identical reference GC outputs across runs.
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Core"))]
    pub seed: Option<u64>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: RefGCWindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Optional BED file(s) with blacklisted regions [path]
    ///
    /// We count no fragment intervals that overlap a blacklisted base.
    /// This results in a lower count for long fragment lengths, which
    /// is not a problem due to length-wise normalization in the downstream
    /// `cfdna gc-bias` command.
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

    /// Whether to skip the smoothing of raw GC counts `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Smoothing"))]
    pub skip_smoothing: bool,

    /// Size of tiles to process the reference in `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "10000000",
            value_parser = clap::value_parser!(u32).range(1000000..),
            help_heading = "Core"
        )
    )]
    pub tile_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl RefGCBiasConfig {
    pub fn check_smoothing_settings(&self) -> anyhow::Result<()> {
        if self.skip_smoothing {
            return Ok(());
        }
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

#[cfg(feature = "cli")]
fn parse_positive_usize(raw: &str) -> Result<usize, String> {
    let value = raw
        .parse::<usize>()
        .map_err(|err| format!("invalid integer: {err}"))?;
    if value == 0 {
        return Err("must be greater than zero".to_string());
    }
    Ok(value)
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct RefGCWindowsArgs {
    /// BED file with regions to include `[path]`
    ///
    /// We count at the **unique positions** included in the specified intervals.
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Windows"))]
    pub by_bed: Option<PathBuf>,
}

impl RefGCWindowsArgs {
    pub fn resolve_windows(&self) -> WindowSpec {
        if let Some(p) = self.by_bed.clone() {
            WindowSpec::Bed(p)
        } else {
            WindowSpec::Global
        }
    }
}
