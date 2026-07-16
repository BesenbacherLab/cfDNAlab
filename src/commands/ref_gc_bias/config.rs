use crate::commands::cli_common::*;
use crate::{ToCliCommand, cli_command::helpers::*};
use std::path::PathBuf;

const DEFAULT_N_POSITIONS: usize = 500_000_000;
const DEFAULT_END_OFFSET: u8 = 10;
const DEFAULT_SMOOTHING_SIGMA: f64 = 0.8;
const DEFAULT_SMOOTHING_RADIUS: u8 = 2;
const DEFAULT_TILE_SIZE: u32 = 10_000_000;

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
#[derive(Debug, Clone, PartialEq)]
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
    ///   `<prefix>.ref_gc_package.zarr`
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Number of threads to use (increases RAM usage) `[integer]`
    ///
    /// Defaults to the number of available CPU cores (-1).
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 't',
            long,
            default_value_t = crate::shared::thread_pool::default_thread_count(),
            help_heading = "Core"
        )
    )]
    pub n_threads: usize,

    /// Number of genomic starting positions to sample `[integer]`
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
            default_value_t = DEFAULT_N_POSITIONS,
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
        clap(long, default_value_t = DEFAULT_END_OFFSET,
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
        clap(long, default_value_t = DEFAULT_SMOOTHING_SIGMA,
             value_parser = clap::value_parser!(f64), help_heading="Smoothing"))]
    pub smoothing_sigma: f64,

    /// Radius of Gaussian kernel that smoothes raw GC counts for each fragment length `[integer]`
    ///
    /// Kernel size is `2 * radius + 1`.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_SMOOTHING_RADIUS,
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
            default_value_t = DEFAULT_TILE_SIZE,
            value_parser = clap::value_parser!(u32).range(1000000..),
            help_heading = "Core"
        )
    )]
    pub tile_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl RefGCBiasConfig {
    /// Build a `ref-gc-bias` config with the same defaults used by the CLI.
    pub fn new(ref_2bit: PathBuf, output_dir: PathBuf, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ref_genome: Ref2BitRequiredArgs { ref_2bit },
            output_dir,
            output_prefix: String::new(),
            n_threads: crate::shared::thread_pool::default_thread_count(),
            n_positions: DEFAULT_N_POSITIONS,
            seed: None,
            windows: RefGCWindowsArgs::default(),
            chromosomes,
            blacklist: None,
            fragment_lengths: FragmentLengthArgs::default(),
            end_offset: DEFAULT_END_OFFSET,
            skip_interpolation: false,
            smoothing_sigma: DEFAULT_SMOOTHING_SIGMA,
            smoothing_radius: DEFAULT_SMOOTHING_RADIUS,
            skip_smoothing: false,
            tile_size: DEFAULT_TILE_SIZE,
            logging: LoggingArgs::default(),
        }
    }

    /// Set the 2bit reference genome path.
    pub fn set_ref_2bit(&mut self, ref_2bit: PathBuf) {
        self.ref_genome.ref_2bit = ref_2bit;
    }

    /// Set the output directory.
    pub fn set_output_dir(&mut self, output_dir: PathBuf) {
        self.output_dir = output_dir;
    }

    /// Set the optional filename prefix for the output package.
    pub fn set_output_prefix<S: Into<String>>(&mut self, output_prefix: S) {
        self.output_prefix = output_prefix.into();
    }

    /// Set the number of worker threads.
    pub fn set_n_threads(&mut self, n_threads: usize) {
        self.n_threads = n_threads;
    }

    /// Set the approximate number of genomic starting positions to sample.
    pub fn set_n_positions(&mut self, n_positions: usize) {
        self.n_positions = n_positions;
    }

    /// Set the optional sampling seed.
    pub fn set_seed(&mut self, seed: Option<u64>) {
        self.seed = seed;
    }

    /// Set the reference-window selection.
    pub fn set_windows(&mut self, windows: RefGCWindowsArgs) {
        self.windows = windows;
    }

    /// Set the optional BED file used to restrict sampled reference positions.
    pub fn set_by_bed(&mut self, by_bed: Option<PathBuf>) {
        self.windows.by_bed = by_bed;
    }

    /// Set the chromosome selection.
    pub fn set_chromosomes(&mut self, chromosomes: ChromosomeArgs) {
        self.chromosomes = chromosomes;
    }

    /// Set optional blacklist BED files.
    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

    /// Return mutable access to the fragment length range.
    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    /// Set the fragment length range.
    pub fn set_fragment_lengths(&mut self, fragment_lengths: FragmentLengthArgs) {
        self.fragment_lengths = fragment_lengths;
    }

    /// Set how many bases to ignore at each fragment end for GC counting.
    pub fn set_end_offset(&mut self, end_offset: u8) {
        self.end_offset = end_offset;
    }

    /// Set whether zero-count interpolation should be skipped.
    pub fn set_skip_interpolation(&mut self, skip_interpolation: bool) {
        self.skip_interpolation = skip_interpolation;
    }

    /// Set the Gaussian smoothing sigma.
    pub fn set_smoothing_sigma(&mut self, smoothing_sigma: f64) {
        self.smoothing_sigma = smoothing_sigma;
    }

    /// Set the Gaussian smoothing radius.
    pub fn set_smoothing_radius(&mut self, smoothing_radius: u8) {
        self.smoothing_radius = smoothing_radius;
    }

    /// Set whether raw GC-count smoothing should be skipped.
    pub fn set_skip_smoothing(&mut self, skip_smoothing: bool) {
        self.skip_smoothing = skip_smoothing;
    }

    /// Set the reference tile size.
    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    /// Set logging options used when rendering equivalent CLI calls.
    pub fn set_logging(&mut self, logging: LoggingArgs) {
        self.logging = logging;
    }

    /// Validate smoothing settings before command execution.
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
#[derive(Debug, Clone, Default, PartialEq)]
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

impl ToCliCommand for RefGCBiasConfig {
    fn to_cli_args(&self) -> crate::Result<Vec<std::ffi::OsString>> {
        let mut args = command_args("ref-gc-bias");
        push_ref_2bit_required(&mut args, &self.ref_genome);
        push_path(&mut args, "--output-dir", &self.output_dir);
        push_output_prefix(&mut args, &self.output_prefix);
        push_value(&mut args, "--n-threads", self.n_threads);
        push_value(&mut args, "--n-positions", self.n_positions);
        if let Some(seed) = self.seed {
            push_value(&mut args, "--seed", seed);
        }
        push_optional_path(&mut args, "--by-bed", self.windows.by_bed.as_deref());
        push_chromosomes(&mut args, &self.chromosomes);
        push_path_values(&mut args, "--blacklist", self.blacklist.as_deref());
        push_fragment_lengths(&mut args, &self.fragment_lengths);
        push_value(&mut args, "--end-offset", self.end_offset);
        push_bool(&mut args, "--skip-interpolation", self.skip_interpolation);
        push_value(&mut args, "--smoothing-sigma", self.smoothing_sigma);
        push_value(&mut args, "--smoothing-radius", self.smoothing_radius);
        push_bool(&mut args, "--skip-smoothing", self.skip_smoothing);
        push_value(&mut args, "--tile-size", self.tile_size);
        push_logging(&mut args, &self.logging);
        Ok(args)
    }
}

#[cfg(test)]
mod tests {
    include!("config_tests.rs");
}
