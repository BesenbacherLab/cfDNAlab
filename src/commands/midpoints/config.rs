use crate::{
    commands::cli_common::{
        ApplyGCArgs, ChromosomeArgs, IOCArgs, LoggingArgs, ScaleGenomeArgs, UnpairedArgs,
        resolve_length_bin_edges,
    },
    shared::{
        blacklist::BlacklistStrategy,
        constants::{MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION},
    },
};
use anyhow::Result;
use std::path::PathBuf;

/// Count positional fragment **midpoint** coverage in groups of genomic windows.
///
/// **Midpoints**: The center of the fragment span, with ties (in even-sized fragments)
/// randomly and reproducibly assigned to the left or right mid-position to avoid bias
/// from always rounding in the same direction.
///
/// **Groups**: The coverage profiles are collapsed (summed per position) across windows in a group.
/// E.g., transcription factors as groups with binding sites as windows, yielding the
/// overall midpoint profile per transcription factor.
///
/// ## Fragment span definition
///
/// **Paired-end**: `[forward.pos, reverse.reference_end)`, the reference span
/// from the first aligned position on the forward read to the last aligned
/// position on the reverse read.
///
/// **Unpaired** where each read is a fragment: `[read.pos, read.reference_end)`,
/// the reference span from the first to the last aligned position on the read.
///
/// The utilized fragment length range is specified via `--length-bins`.
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
pub struct MidpointsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls of the command.
    ///
    /// Examples produce files like:
    ///   `<prefix>.midpoint_profiles.npy`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// The grouped fixed-size intervals to count within `[path]`
    ///
    /// A BED-like file of genomic intervals and their respective group names.
    ///
    /// Must be sorted by the `chromosome` and `start` coordinates, and
    /// all intervals must have the same size.
    ///
    /// Sites with the same group name are collapsed to a single profile.
    ///
    /// Columns: `chromosome, start, end, group_name`. No header.
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'w',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub intervals: PathBuf,

    /// Edges of fragment length bins to count in `[string(s)]`
    ///
    /// This also defines the min and max fragment lengths.
    ///
    /// Accepted forms:
    ///
    /// - A single value with `start:end:step`:
    ///   Creates contiguous bins from `start` to `end` (end-exclusive) in `step` increments.
    ///   Example: `30:1000:10` -> bins `[30,40), [40,50), ..., [990,1000)`.
    ///
    /// - Multiple integer values interpreted as bin edges:
    ///   Example: `--length-bins 30 80 150 220 500 1001` -> bins `[30,80), [80,150), ..., [500,1001)`.
    ///
    /// **NOTE**: Memory consumption increases linearly with the number of bins.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            default_values_t = [String::from("30"), String::from("1001")],
            help_heading = "Core"
        )
    )]
    pub length_bins: Vec<String>,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    ///
    /// This is **NOT** recommended by default as it trims the tails of the length distribution.
    ///
    /// Note, that we only keep inward-directed fragments within the specified length range, so
    /// there's no real need for proper-pair filtering.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// **NOTE**: It may be an advantage to instead remove intervals that lie within
    /// half the maximum fragment length of blacklisted regions from the `--intervals` file.
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
    /// `midpoint` checks the single central base for odd fragments and either
    /// central base for even fragments.
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

    /// Group indices to plot as midpoint profiles `[integers]`
    ///
    /// Comma separated list of zero-based group indices to plot after counting.
    ///
    /// This plotting step is intended for quick QC of the outputs. It's not
    /// optimized for publication etc. (although feel free!)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_delimiter = ',',
            default_values_t = [0_usize],
            help_heading = "Plotting"
        )
    )]
    pub plot_groups: Vec<usize>,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl MidpointsConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs, intervals: PathBuf) -> Self {
        Self {
            ioc,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix: String::new(),
            intervals,
            length_bins: vec!["30".to_string(), "1001".to_string()],
            tile_size: 20000000,
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
            gc: ApplyGCArgs {
                gc_file: None,
                gc_tag: None,
                neutralize_invalid_gc: false,
            },
            ref_2bit: None,
            plot_groups: vec![0],
            logging: LoggingArgs::default(),
        }
    }

    pub fn set_output_prefix<S: Into<String>>(&mut self, prefix: S) {
        self.output_prefix = prefix.into();
    }

    pub fn set_length_bins(&mut self, edges: Vec<u32>) {
        assert!(
            edges.len() >= 2,
            "length bin edges must contain at least two values"
        );
        self.length_bins = edges.into_iter().map(|edge| edge.to_string()).collect();
    }

    pub fn set_length_bins_spec<S: Into<String>>(&mut self, spec: S) {
        self.length_bins = vec![spec.into()];
    }

    pub fn resolve_length_bins(&self) -> Result<Vec<u32>> {
        resolve_length_bin_edges(
            &self.length_bins,
            MIN_ACGT_BASES_FOR_GC_FRACTION,
            MAX_SUPPORTED_FRAGMENT_LENGTH,
        )
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.scale_genome = scale;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}

#[cfg(test)]
mod tests {
    include!("config_tests.rs");
}
