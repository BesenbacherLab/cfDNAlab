use crate::commands::cli_common::{
    ApplyGCArgFileOnly, ChromosomeArgs, FragmentLengthArgs, ScaleGenomeArgs, UnpairedArgs,
    WindowSpec,
};
use crate::shared::blacklist::BlacklistStrategy;
use std::path::PathBuf;

/// Apply filtering and corrections to the fragments in a BAM file
/// and write to a new coordinate-sorted BAM file.
///
/// To use our corrections and filters in *custom*, downstream analyses, you can apply
/// them directly to a given BAM file. Filter which reads/fragments to write and add correction
/// weights as AUX tags on the reads. The new BAM file is coordinate-sorted.
///
/// **NOTE**: This is **not** needed for running other `cfDNAlab` tools.
/// Those tools will **not** automatically use the correction tags.
///
/// ## Genomic smoothing (--scaling-factors)
///
/// The coverage weight that would normally be **multiplied** with the fragment's count value (`1.0`)
/// is written as the AUX tag `COV` in the read(s).
///
/// ## GC bias correction
///
/// The GC bias correction weight that would normally be **multiplied** with the fragment's count
/// value (`1.0` or the smoothed value) is written as the AUX tag `GC` in the read(s).
///
/// ## Fragment length
///
/// The fragment length is written to the AUX tag "FLEN".
///
/// Definition:
///
/// **Paired-end**: `end(reverse) - start(forward)`.
///
/// **Unpaired** where each read is a fragment: `end(read) - start(read)`.
///
/// ## Always-on exclusion criteria
///
/// The following criteria always exclude a read:
///
/// The read is secondary, supplementary or duplicate.
/// The read failed quality check.
///
/// **Paired-end input only**:
/// The read or mate read is unmapped.
/// The read is mapped to a different `tid` than the mate.
/// The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct BamToBamConfig {
    /// Indexed, coordinate-sorted BAM input file `[path]`
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
    pub in_bam: PathBuf,

    /// Path to write coordinate-sorted BAM at `[path]`
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
    pub out_bam: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Intervals to keep overlapping fragments from `[path]`
    ///
    /// Reads that are part of a fragment that overlaps a window
    /// are considered for the new BAM file.
    #[cfg_attr(
        feature = "cli",
        clap(long = "by-bed", value_parser, help_heading = "Windows")
    )]
    pub by_bed: Option<PathBuf>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Keep the specified chromosome order instead of sorting lexicographically `[flag]`
    ///
    /// Many tools expect BAM files to be sorted as `chr1, chr10, chr11, ...`. By default,
    /// we thus sort the specified chromosomes lexicographically. This is different to other
    /// commands in `cfDNAlab`, which directly use the passed order of chromosomes.
    #[cfg_attr(
        feature = "cli",
        clap(long = "skip-chromosome-sort", help_heading = "Core")
    )]
    pub skip_chromosome_sort: bool,

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
    pub gc: ApplyGCArgFileOnly,

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

impl BamToBamConfig {
    pub fn new(in_bam: PathBuf, out_bam: PathBuf, chromosomes: ChromosomeArgs) -> Self {
        Self {
            in_bam,
            out_bam,
            by_bed: None,
            chromosomes,
            skip_chromosome_sort: false,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs::default(),
            min_mapq: 0,
            require_proper_pair: false,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::Any,
            gc: ApplyGCArgFileOnly {
                gc_file: None,
                drop_invalid_gc: false,
            },
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

    pub fn set_by_bed(&mut self, by_bed: Option<PathBuf>) {
        self.by_bed = by_bed;
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.scale_genome = scale;
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
    pub fn set_skip_chromosome_sort(&mut self, skip_chromosome_sort: bool) {
        self.skip_chromosome_sort = skip_chromosome_sort;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgFileOnly) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}
