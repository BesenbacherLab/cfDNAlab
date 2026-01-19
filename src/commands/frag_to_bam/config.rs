use crate::commands::cli_common::{ChromosomeArgs, FragmentLengthArgs};
use crate::shared::blacklist::BlacklistStrategy;
use std::path::PathBuf;

/// Convert a finaleDB-style frag file to a BAM file with unpaired reads
/// (each read is a full fragment).
///
/// The first five columns in the frag file:
/// `Chromosome, Start, End, MapQ, Strand`.
///
/// Other columns are ignored.
///
/// Each read in the new BAM file represents a fragment from the frag file.
///
/// The BAM header contains all contigs from `--chrom-sizes` in the `--chrom-sizes` order.
///
/// The BAM file is not indexed. This can be done with `samtools index`.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
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

    /// Prefix for output file (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.bam`,
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value = "fragments", help_heading = "Core")
    )]
    pub output_prefix: String,

    /// File with chromosome sizes (FAI or two-column sizes) for the BAM header `[path]`
    ///
    /// E.g. the USCS `hg38.chrom.sizes` file (or similar for your assembly).
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
            output_prefix: "fragments".into(),
            chromosomes,
            chrom_sizes,
            fragment_lengths: FragmentLengthArgs::default(),
            min_mapq: 0,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::Any,
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

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
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
