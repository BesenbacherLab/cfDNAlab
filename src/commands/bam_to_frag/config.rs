use crate::commands::cli_common::{
    ApplyGCArgs, ChromosomeArgs, FragmentLengthArgs, IOCArgs, ScaleGenomeArgs, WindowSpec,
};
use crate::shared::blacklist::BlacklistStrategy;
use std::path::PathBuf;

/// Write the fragments from a BAM file to a finaleDB-style frag file.
///
/// Information in the `.frag.tsv` file:
///
///  - **Chromosome**
///
///  - **Start**: forward.pos
///
///  - **End**: reverse.end
///
///  - **MapQ**: Minimum mapping quality for the two reads
///
///  - **Strand**: The strand alignment of read1
///
/// AND, when either/both `--gc-file` and `--scaling-factors` are specified:
///
///  - **GC Weight**: The multiplicative weight needed to correct for GC bias.
///
///  - **Scaling Weight**: The multiplicative weight needed to perform genomic smoothing.
///
/// Note: When GC correction is not specified but genomic scaling is, the sixth column is the scaling weight.
///
/// The accompanying `*.frag.header.tsv` file has the matching column names.
///
/// Fragments are sorted by `(chromosome, start, end)`, using the chromosome order in `--chromosomes`.
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
#[derive(Clone)]
pub struct BamToFragConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    /// Prefix for output file (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.frag.tsv.gz`,
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value = "fragments", help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Intervals to keep overlapping fragments from `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(long = "by-bed", value_parser, help_heading = "Windows")
    )]
    pub by_bed: Option<PathBuf>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include `[integer]`
    ///
    /// Defaults to 0 to allow making filtering decisions downstream.
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "0", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    ///
    /// This is NOT recommended by default, as it trims the tails of the length distribution.
    /// It may be useful to match the files in FinaleDB.
    ///
    /// Note, that we only keep inward-directed fragments within a specified length range, so
    /// there's no real need for proper-pair filtering.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

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
    pub gc: ApplyGCArgs,

    /// Optional 2bit reference genome file [path]
    ///
    /// NOTE: Required for GC correction, otherwise ignored.
    ///
    /// E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'r',
            long,
            value_parser,
            required = false,
            help_heading = "GC Correction"
        )
    )]
    pub ref_2bit: Option<PathBuf>,
}

impl BamToFragConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            output_prefix: "coverage".into(),
            by_bed: None,
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs::default(),
            min_mapq: 0,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::Any,
            gc: ApplyGCArgs { gc_file: None },
            ref_2bit: None,
        }
    }

    /// If neither flag is set, default to `Global`.
    pub fn resolve_windows(&self) -> WindowSpec {
        if let Some(p) = self.by_bed.clone() {
            WindowSpec::Bed(p)
        } else {
            WindowSpec::Global
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_by_bed(&mut self, by_bed: Option<PathBuf>) {
        self.by_bed = by_bed;
    }

    pub fn set_scale_genome(&mut self, scale_genome: ScaleGenomeArgs) {
        self.scale_genome = scale_genome;
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
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

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}
