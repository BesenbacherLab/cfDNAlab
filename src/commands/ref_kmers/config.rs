use crate::commands::cli_common::*;
use crate::{ToCliCommand, cli_command::helpers::*};
use anyhow::{Result, ensure};
use std::path::PathBuf;

const DEFAULT_TILE_SIZE: u32 = 10_000_000;

/// Count reference k-mer frequencies for genomic windows or groups.
///
/// Builds a reference-sequence background for downstream k-mer correction.
/// It writes row-wise frequencies plus a row scaling factor
/// for reconstructing counts downstream.
///
/// K-mers are counted left-to-right and only contains the forward-oriented motifs
/// from eligible positions (non-blacklisted regions).
/// Downstream applications needing the reverse-oriented motifs would need to
/// reverse-complement the motifs.
///
/// The selected window assignment controls whether overlapping k-mers contribute
/// fractional counts (default) or full counts when overlapping with a certain proportion.
///
/// K-mers overlapping blacklisted bases or ambiguous `N` bases are not counted.
///
/// ## Limitations
///
/// For **end-motifs** (see `cfdna ends`), beware that this command also includes
/// k-mers at the edges of chromosomes and blacklisted regions, where one of the fragment's
/// ends is not actually countable. The counts do thus not fully represent the "possible
/// fragment ends", although the error is tiny in non-small windows compared
/// to a per-fragment-length possible-fragments count.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, PartialEq)]
pub struct RefKmersConfig {
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
    ///   `<prefix>.ref_kmer_counts.zarr`
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

    /// Size of the k-mers `[integer]`
    ///
    /// Without `--motifs-file`, the largest supported k-mer size is `27`.
    ///
    /// With `--motifs-file`, larger k-mers can be counted because the command only tracks the
    /// selected k-mer subspace.
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'k',
            long,
            value_parser = clap::value_parser!(u8).range(1..),
            required = true,
            help_heading = "Core"
        )
    )]
    pub kmer_size: u8,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: DistributionWindowsArgs,

    /// How to assign k-mers to windows `[string]`
    ///
    /// Possible values:
    ///     `"count-overlap"`, `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
    ///
    /// `"count-overlap"`: Count the fraction of k-mer bases overlapping each window.
    ///
    /// `"any"`: Count the k-mer in every window overlapping at least one k-mer base.
    ///
    /// `"all"`: Count the k-mer only in windows containing all k-mer bases.
    ///
    /// `"midpoint"`: Count the k-mer in the window overlapping its center base.
    /// This mode requires an odd `--kmer-size`.
    ///
    /// `"proportion=<threshold>"`: Count the k-mer when at least this fraction of its bases
    /// overlaps a window.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "count-overlap",
            ignore_case = true,
            help_heading = "Window Assignment"
        )
    )]
    pub assign_by: WindowAssigner,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Optional BED file(s) with blacklisted regions [path]
    ///
    /// K-mers overlapping a blacklisted base are not counted.
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'b',
            long,
            value_parser,
            num_args = 1..,
            action = clap::ArgAction::Append,
            help_heading="Filtering"
        )
    )]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Collapse each k-mer with its reverse complement `[flag]`
    ///
    /// Odd-sized k-mers are collapsed such that the middle base is `A` or `C`.
    /// Even-sized k-mers are collapsed to the lexicographically lowest motif.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub canonical: bool,

    /// Include every possible motif in the output, even if its count is zero `[flag]`
    ///
    /// **NOTE**: When `--motifs-file` is specified, it defines the "possible" motifs.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Motifs"))]
    pub all_motifs: bool,

    /// File with motifs to include `[path]`
    ///
    /// TSV-like file (tab-separated, no header) with one motif per line.
    /// Add a second column with a group name to count multiple motifs together.
    ///
    /// Each motif must be an A/C/G/T k-mer of length `--kmer-size`.
    /// For compatibility with `cfdna ends`, a single `_` separator is accepted and removed,
    /// so `AC_GT` is read as `ACGT`.
    ///
    /// In grouped mode, the output motif axis contains one entry per distinct group name,
    /// ordered alphabetically by group name.
    ///
    /// Frequencies are normalized over the selected motifs or groups in this file.
    /// Unlisted k-mers are not part of the denominator. The row scaling factor stores the
    /// selected-target count total, so reconstructed counts are selected k-mer or group counts.
    ///
    /// Specifying the allowed subset of motifs beforehand enables counting of much
    /// larger k-mers without enumerating the full k-mer universe.
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Motifs"))]
    pub motifs_file: Option<PathBuf>,

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

impl RefKmersConfig {
    /// Validate settings that can be checked before reading inputs.
    pub(crate) fn validate(&self) -> Result<()> {
        validate_output_prefix(&self.output_prefix)?;

        let window_opt = self.windows.resolve_windows();
        if let DistributionWindowSpec::Size(window_bp) = &window_opt {
            ensure!(*window_bp > 0, "`--by-size` must be greater than 0");
        }
        ensure!(self.kmer_size > 0, "`--kmer-size` must be greater than 0");
        ensure!(
            !matches!(self.assign_by, WindowAssigner::Midpoint) || self.kmer_size % 2 == 1,
            "`--assign-by midpoint` requires an odd `--kmer-size`"
        );
        ensure!(
            self.motifs_file.is_some()
                || usize::from(self.kmer_size)
                    <= crate::shared::kmers::kmer_codec::MAX_RADIX5_KMER_SIZE,
            "`--kmer-size` > {} requires `--motifs-file`",
            crate::shared::kmers::kmer_codec::MAX_RADIX5_KMER_SIZE
        );

        Ok(())
    }

    /// Build a `ref-kmers` config with the same defaults used by the CLI.
    pub fn new(
        ref_2bit: PathBuf,
        output_dir: PathBuf,
        kmer_size: u8,
        chromosomes: ChromosomeArgs,
    ) -> Self {
        Self {
            ref_genome: Ref2BitRequiredArgs { ref_2bit },
            output_dir,
            output_prefix: String::new(),
            n_threads: crate::shared::thread_pool::default_thread_count(),
            kmer_size,
            windows: DistributionWindowsArgs::default(),
            assign_by: WindowAssigner::default(),
            chromosomes,
            blacklist: None,
            canonical: false,
            all_motifs: false,
            motifs_file: None,
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

    /// Set the k-mer size.
    pub fn set_kmer_size(&mut self, kmer_size: u8) {
        self.kmer_size = kmer_size;
    }

    /// Set the reference-window selection.
    pub fn set_windows(&mut self, windows: DistributionWindowsArgs) {
        self.windows = windows;
    }

    /// Set how k-mers are assigned to windows.
    pub fn set_assign_by(&mut self, assign_by: WindowAssigner) {
        self.assign_by = assign_by;
    }

    /// Set the chromosome selection.
    pub fn set_chromosomes(&mut self, chromosomes: ChromosomeArgs) {
        self.chromosomes = chromosomes;
    }

    /// Set optional blacklist BED files.
    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

    /// Set whether to collapse each k-mer with its reverse complement.
    pub fn set_canonical(&mut self, canonical: bool) {
        self.canonical = canonical;
    }

    /// Set whether zero-count motifs should be included in the output.
    pub fn set_all_motifs(&mut self, all_motifs: bool) {
        self.all_motifs = all_motifs;
    }

    /// Set the optional motif list file.
    pub fn set_motifs_file(&mut self, motifs_file: Option<PathBuf>) {
        self.motifs_file = motifs_file;
    }

    /// Set the reference tile size.
    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    /// Set logging options used when rendering equivalent CLI calls.
    pub fn set_logging(&mut self, logging: LoggingArgs) {
        self.logging = logging;
    }
}

impl ToCliCommand for RefKmersConfig {
    fn to_cli_args(&self) -> crate::Result<Vec<std::ffi::OsString>> {
        let mut args = command_args("ref-kmers");
        push_ref_2bit_required(&mut args, &self.ref_genome);
        push_path(&mut args, "--output-dir", &self.output_dir);
        push_output_prefix(&mut args, &self.output_prefix);
        push_value(&mut args, "--n-threads", self.n_threads);
        push_value(&mut args, "--kmer-size", self.kmer_size);
        push_distribution_windows(&mut args, &self.windows);
        push_value(
            &mut args,
            "--assign-by",
            window_assigner_value(&self.assign_by),
        );
        push_chromosomes(&mut args, &self.chromosomes);
        push_path_values(&mut args, "--blacklist", self.blacklist.as_deref());
        push_bool(&mut args, "--canonical", self.canonical);
        push_bool(&mut args, "--all-motifs", self.all_motifs);
        push_optional_path(&mut args, "--motifs-file", self.motifs_file.as_deref());
        push_value(&mut args, "--tile-size", self.tile_size);
        push_logging(&mut args, &self.logging);
        Ok(args)
    }
}
