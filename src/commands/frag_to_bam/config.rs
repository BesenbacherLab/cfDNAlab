use crate::commands::cli_common::{ChromosomeArgs, FragmentLengthArgs, LoggingArgs};
use crate::shared::blacklist::BlacklistStrategy;
use crate::{ToCliCommand, cli_command::helpers::*};
use std::path::PathBuf;

/// Convert a finaleDB-style frag file to a BAM file with unpaired reads
/// (each read is a full fragment).
///
/// Each read in the new BAM file represents a fragment from the frag file.
///
/// The first five columns in the frag file should be:
/// `Chromosome, Start, End, MapQ, Strand`.
///
/// ## Extra columns
///
/// Optional extra columns can be transferred to BAM AUX tags when column names are known.
///
/// The recognized names and respective AUX tags are:
///
/// - `gc_weight` -> `GC`
///
/// - `coverage_scaling_weight` -> `cw`
///
/// - `count_scaling_weight` -> `nw`
///
/// - `flen` -> `fl`
///
/// ## BAM file
///
/// The BAM header contains all contigs from `--chrom-sizes` in the `--chrom-sizes` order.
///
/// A BAI index is written next to the BAM as `<output>.bam.bai`.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, PartialEq)]
pub struct FragToBamConfig {
    /// Path to a coordinate-sorted `.tsv` frag file `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'i',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub frag: PathBuf,

    /// Output directory to write new BAM file in `[path]`
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

    /// Optional prefix for output file (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.fragments.bam`
    ///
    /// With an empty prefix, the output filename is `fragments.bam`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Optional header file with tab-separated column names for the frag file [path]
    ///
    /// Supply this when you want to transfer extra columns
    /// (`gc_weight`, `coverage_scaling_weight`, `count_scaling_weight`, and/or `flen`) to AUX tags
    /// in the BAM file and the frag file has no inline header row.
    ///
    /// **Auto-detection**: The command also tries to auto-detect a companion header file
    /// named `<prefix>.frag.header.tsv` when the frag path follows
    /// `<prefix>.frag.tsv` (optionally with `.gz` or `.zst`).
    ///
    /// When no headers are supplied/detected or found inline, the command still accepts
    /// headerless 5-column frag files.
    ///
    /// Use `--ignore-extras` when you want to ignore all extra columns after the first five.
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Core"))]
    pub frag_header: Option<PathBuf>,

    /// File with chromosome sizes (FAI or two-column sizes) for the BAM header `[path]`
    ///
    /// E.g. the UCSC `hg38.chrom.sizes` file (or similar for your assembly).
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Core"))]
    pub chrom_sizes: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include `[integer]`
    ///
    /// Defaults to 0 to allow making filtering decisions downstream.
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "0", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Ignore all frag columns after the first five `[flag]`
    ///
    /// This disables mapping extra columns to BAM AUX tags.
    /// It also allows headers with extra names that are not supported for AUX mapping.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub ignore_extras: bool,

    /// Allow unknown extra header columns and ignore them `[flag]`
    ///
    /// By default, unknown extra columns cause an error to prevent silent mistakes.
    ///
    /// With this flag, unknown extra columns are ignored with a warning, while known
    /// extra columns (`gc_weight`, `coverage_scaling_weight`, `count_scaling_weight`, `flen`) are still transferred.
    ///
    /// If you want to ignore all extras, use `--ignore-extras` instead.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub allow_unknown_extras: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-min-size",
            default_value = "1",
            help_heading = "Filtering"
        )
    )]
    pub blacklist_min_size: u64,

    /// The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
    ///
    /// Possible values:
    ///     `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
    ///
    /// Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-strategy",
            default_value = "any",
            ignore_case = true,
            help_heading = "Filtering"
        )
    )]
    pub blacklist_strategy: BlacklistStrategy,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl FragToBamConfig {
    pub fn new(
        frag: PathBuf,
        output_dir: PathBuf,
        chromosomes: ChromosomeArgs,
        chrom_sizes: PathBuf,
    ) -> Self {
        Self {
            frag,
            output_dir,
            output_prefix: String::new(),
            chromosomes,
            chrom_sizes,
            frag_header: None,
            fragment_lengths: FragmentLengthArgs::default(),
            min_mapq: 0,
            ignore_extras: false,
            allow_unknown_extras: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::Any,
            logging: LoggingArgs::default(),
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_frag(&mut self, frag: PathBuf) {
        self.frag = frag;
    }

    pub fn set_output_dir(&mut self, output_dir: PathBuf) {
        self.output_dir = output_dir;
    }

    pub fn set_chromosomes(&mut self, chromosomes: ChromosomeArgs) {
        self.chromosomes = chromosomes;
    }

    pub fn set_chrom_sizes(&mut self, chrom_sizes: PathBuf) {
        self.chrom_sizes = chrom_sizes;
    }

    pub fn set_frag_header(&mut self, frag_header: Option<PathBuf>) {
        self.frag_header = frag_header;
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_ignore_extras(&mut self, ignore_extras: bool) {
        self.ignore_extras = ignore_extras;
    }

    pub fn set_allow_unknown_extras(&mut self, allow_unknown_extras: bool) {
        self.allow_unknown_extras = allow_unknown_extras;
    }

    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

    pub fn set_blacklist_min_size(&mut self, blacklist_min_size: u64) {
        self.blacklist_min_size = blacklist_min_size;
    }

    pub fn set_blacklist_strategy(&mut self, blacklist_strategy: BlacklistStrategy) {
        self.blacklist_strategy = blacklist_strategy;
    }
}

impl ToCliCommand for FragToBamConfig {
    fn to_cli_args(&self) -> crate::Result<Vec<std::ffi::OsString>> {
        let mut args = command_args("frag-to-bam");
        push_path(&mut args, "--frag", &self.frag);
        push_path(&mut args, "--output-dir", &self.output_dir);
        push_output_prefix(&mut args, &self.output_prefix);
        push_optional_path(&mut args, "--frag-header", self.frag_header.as_deref());
        push_path(&mut args, "--chrom-sizes", &self.chrom_sizes);
        push_chromosomes(&mut args, &self.chromosomes);
        push_fragment_lengths(&mut args, &self.fragment_lengths);
        push_value(&mut args, "--min-mapq", self.min_mapq);
        push_bool(&mut args, "--ignore-extras", self.ignore_extras);
        push_bool(
            &mut args,
            "--allow-unknown-extras",
            self.allow_unknown_extras,
        );
        push_blacklist_common(
            &mut args,
            self.blacklist.as_deref(),
            self.blacklist_min_size,
            &self.blacklist_strategy,
        );
        push_logging(&mut args, &self.logging);
        Ok(args)
    }
}
