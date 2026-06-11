use crate::{
    commands::{
        cli_common::{
            ApplyGCArgFileOnly, AssignToWindowArgs, ChromosomeArgs, DistributionWindowsArgs,
            IOCArgs, LoggingArgs, ScaleGenomeArgs, UnpairedArgs, resolve_length_bin_edges,
        },
        gc_bias::correct::{GCLengthRange, MarginalizeLengthsWeightingScheme},
    },
    shared::{
        blacklist::BlacklistStrategy,
        clip_mode::ClipMode,
        constants::{
            DEFAULT_MAX_SOFT_CLIPS, MAX_MAX_SOFT_CLIPS, MAX_SUPPORTED_FRAGMENT_LENGTH,
            MIN_ACGT_BASES_FOR_GC_FRACTION,
        },
        indel_mode::IndelMode,
    },
};
use anyhow::Result;
use std::path::PathBuf;

pub const DEFAULT_MAX_DELETION_BASES: u16 = 100;
pub const MAX_DELETION_BASES: u16 = 256;
pub const DEFAULT_OUTPUT_DECIMALS: u8 = 6;

/// Count fragment lengths in a BAM-file.
///
/// Writes a wide compressed TSV count table. Rows contain the global output,
/// genomic windows, or grouped-BED groups depending on the selected windowing
/// mode. Single-bp fragment length bins are stored as `count_<length>` columns.
/// Wider bins are stored as `count_<start>_<end>` columns, where `<start>` and
/// `<end>` are the half-open fragment length bin bounds.
///
/// ## Fragment length definition
///
/// **Paired-end**: `end(reverse) - start(forward)`.
///
/// **Unpaired** where each read is a fragment: `end(read) - start(read)`.
///
/// See also `--indel-mode` and `--clip-mode` for adjusting the length to the
/// present indels and soft clips. When enabled, fragment length filtering is based
/// on the adjusted length.
///
/// ## GC correction
///
/// Weight the contribution of each fragment based on their GC contents.
///
/// The length-dimension of the original correction matrix is averaged out over
/// `--gc-length-range` with a specifiable weighting scheme (see `--gc-length-weighting`).
///
/// The GC percentage is calculated from the **aligned** reference span.
/// It does not consider `--indel-mode` or `--clip-mode`.
///
/// ## Genomic smoothing (--scaling-factors)
///
/// Weight how genomic regions contribute to the length distribution(s), e.g., to reduce the
/// influence of copy number alterations. This weights the contribution of each fragment
/// by region-wise precomputed scaling factors.
///
/// Can be precomputed with `cfdna fragment-count-weights` (recommended) or `cfdna coverage-weights`.
///
/// ## Window assignment
///
/// By default, fragments are counted by their window-overlap fraction. That is, most
/// fragments are counted as `1.0` (before correction/scaling), while fragments overlapping the
/// edge of a window are counted as the fraction it overlaps the window (`< 1.0`).
///  
/// For consecutive non-overlapping windows, this conserves the total mass, as an edge-overlapping
/// fragment will count `f` in one window and `1-f` in the other window.
///
/// With the default width-1 length bins, you can convert counts to base-weighted counts
/// (i.e., the coverage in the window) by multiplying each column by the fragment length
/// it represents. Remember to account for the minimum fragment length offset.
///
/// Other options include counting the full fragment if the **fragment midpoint** or a given
/// **proportion** of positions overlaps the window.
///
/// ## Blacklisting
///
/// Ignores fragments that overlap blacklisted regions with a given proportion.
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
pub struct LengthsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.length_counts.tsv.zst`
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Decimals to round count values to when writing `[integer]`
    ///
    /// This only affects the text representation in the final output.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value_t = DEFAULT_OUTPUT_DECIMALS, value_parser = clap::value_parser!(u8).range(0..), help_heading="Core"))]
    pub decimals: u8,

    /// How to handle insertions and deletions in fragments `[string]`
    ///
    /// Deletions: Both 'D' and 'N' in the cigar string are considered deletions.
    ///
    /// Possible values:
    ///
    /// - `"ignore"`:
    ///   Ignore whether indels are present or not.
    ///
    ///   Lengths are calculated from the reference coordinates `end(reverse) - start(forward)`.
    ///
    /// - `"adjust"`:
    ///   Adjust the reference length by the observed insertions and deletions
    ///   (we cannot adjust in the mate-gap).
    ///
    ///   For bases only covered by a single read, all insertions and deletions are adjusted for.
    ///
    ///   In the mate-overlap, only adjust when both reads show the indel at the same reference position.
    ///   
    ///   Deletions: subtract the reference bases deleted in both reads.
    ///     
    ///   Insertions: add the shortest insertion length per position.
    ///
    ///   **NOTE**: Blacklist exclusion, GC correction, and calculation of scaling weights
    ///   (--scaling-factors) use the aligned reference span.
    ///
    /// - `"skip"`:
    ///   Skip fragments with any insertion or deletion present.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "ignore",
            ignore_case = true,
            help_heading = "Indels and clipping"
        )
    )]
    pub indel_mode: IndelMode,

    /// Skip fragments with more deleted reference bases than this
    /// **when using** `--indel-mode adjust` `[integer]`
    ///
    /// Both `D` and `N` CIGAR operations count as deletion bases.
    ///
    /// **NOTE**: This cap is only used with `--indel-mode adjust`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_MAX_DELETION_BASES,
            value_parser = clap::value_parser!(u16).range(0..=MAX_DELETION_BASES as i64),
            help_heading = "Indels and clipping"
        )
    )]
    pub max_deletion_bases: u16,

    /// How to handle soft clipping in fragment ends `[string]`
    ///
    /// When you believe soft clipping in the fragment ends is mostly due
    /// to alignment difficulties instead of technical artefacts, you can
    /// include the clipped bases in the fragment length.
    ///
    /// Possible values:
    ///
    /// - `"aligned"`:
    ///   Ignore clipped bases and use the aligned positions.
    ///
    /// - `"adjust"`:
    ///   Adjust the fragment length by the observed soft clipped bases in the fragment ends.
    ///
    ///   For paired-end data, the clipping is only considered
    ///   for the 5' ends (start(forward), end(reverse)).
    ///
    ///   **NOTE**: Blacklist exclusion, GC correction, and scaling weights
    ///   (--scaling-factors) use the aligned reference span.
    ///   When `--assign-by count-overlap`, clipped-only window contributions use
    ///   the nearest aligned reference base for scaling.
    ///   
    ///
    /// - `"skip"`:
    ///   Skip fragments with any clipping.
    ///
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "aligned",
            ignore_case = true,
            help_heading = "Indels and clipping"
        )
    )]
    pub clip_mode: ClipMode,

    /// Skip fragments where one or both ends have more soft-clipped bases than
    /// this **when using** `--clip-mode adjust` `[integer]`
    ///
    /// Use `--clip-mode skip` to discard all soft-clipped fragments.
    ///
    /// **NOTE**: This cap is only used with `--clip-mode adjust`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_MAX_SOFT_CLIPS,
            value_parser = clap::value_parser!(u16).range(0..=MAX_MAX_SOFT_CLIPS as i64),
            help_heading = "Indels and clipping"
        )
    )]
    pub max_soft_clips: u16,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: DistributionWindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub window_assignment: AssignToWindowArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    /// Edges of fragment length bins to count in `[string(s)]`
    ///
    /// This also defines the minimum and maximum included fragment lengths.
    ///
    /// Bins are half-open. For example, `--length-bins 10 151 221` creates
    /// bins `[10,151)` and `[151,221)`.
    ///
    /// Accepted forms:
    ///
    /// - A single value with `start:end:step`:
    ///   Creates contiguous bins from `start` to `end` in `step` increments.
    ///   Example: The default `30:1001:1` creates one column per length from 30 through 1000.
    ///
    /// - Multiple integer values interpreted as bin edges:
    ///   Example: `--length-bins 30 80 151 221` creates bins `[30,80)`,
    ///   `[80,151)`, and `[151,221)`.
    ///
    /// **NOTE**: Memory consumption increases linearly with the number of bins.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            default_values_t = [String::from("30:1001:1")],
            help_heading = "Core"
        )
    )]
    pub length_bins: Vec<String>,

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

    /// How to weight the GC-package fragment length bins when estimating the
    /// length-agnostic GC correction `[string]`
    ///
    /// The default GC correction package stores a `fragment length x GC` matrix with
    /// one normalized GC correction curve per fragment length bin. **NOTE**: These
    /// fragment length bins are independent of those specified on `--length-bins`
    /// and are referred to as **GC length bins**.
    ///
    /// When estimating the fragment length distribution itself, using those normalized
    /// correction curves directly would just preserve the original length distribution
    /// (when applied to the same fragments seen by `cfdna gc-bias`).
    ///
    /// We therefore average out the fragment length dimension to get a single, length-agnostic GC correction curve.
    ///
    /// First, `--gc-length-range` selects which GC length bins to use. Then,
    /// `--gc-length-trim-rare` can exclude a fraction of the selected bins
    /// with the lowest frequencies in the length distribution used to build
    /// the GC correction package. Finally, `--gc-length-weighting` controls
    /// how the retained correction curves are collapsed.
    ///
    /// We have three weighting options:
    ///
    /// - `"equal"` weighting (default): Weight every GC length bin the same.
    ///
    ///   Keeps the correction independent of the length distribution we are trying to estimate.
    ///
    ///   Downside: Rare GC length bins contribute the same as the most common GC length bins.
    ///   
    ///   For low-coverage BAM files, this could make the correction more volatile to outliers.
    ///   To reduce this effect, `--gc-length-trim-rare` allows filtering out a fraction of the rarest bins.
    ///
    /// - `"frequency"` weighting: Weight GC length bins by their frequencies in the length distribution
    ///   used to build the GC package.
    ///
    ///   This makes the collapsed curve represent the fragments most often seen in the package.
    ///
    ///   Downside: **Biases** the correction based on the length distribution we are trying to estimate.
    ///   This is **circular**: the length signal is partly used to correct itself.
    ///
    /// - `"max-frequency"` weighting: Use the GC correction curve for the GC length bin with highest
    ///   frequency in the length distribution used to build the GC package.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "equal",
            ignore_case = true,
            help_heading = "GC Correction"
        )
    )]
    pub gc_length_weighting: MarginalizeLengthsWeightingScheme,

    /// Which GC-package fragment length bins to use when averaging out the GC correction matrix `[string]`
    ///
    /// The GC correction package stores one GC correction curve per fragment length bin
    /// (separate from `--length-bins` and referred to as **GC length bins**).
    /// This option controls which GC length bins `--gc-length-weighting`
    /// collapses to a single length-agnostic GC curve.
    ///
    /// Possible values:
    ///
    /// - `"requested"`: Use GC length bins that overlap the range requested by
    ///   `--length-bins`.
    ///
    /// - `"package"`: Use all GC length bins in the GC correction package.
    ///   Useful for repeated calls with different length ranges that you wish to compare downstream.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "requested",
            ignore_case = true,
            help_heading = "GC Correction"
        )
    )]
    pub gc_length_range: GCLengthRange,

    /// Exclude low-frequency selected GC length bins up to this cumulative fraction `[float]`
    ///
    /// After selecting the GC length bins with `--gc-length-range`,
    /// this excludes the least frequent GC length bins while keeping
    /// at least a `1 - fraction` total frequency fraction.
    ///
    /// Conservative: GC length bins with practically the same frequency are grouped,
    /// and the whole group is excluded only if the retained sum of frequencies
    /// is at least `1 - fraction` of the selected total.
    ///
    /// Use this with `--gc-length-weighting equal` when almost-unobserved
    /// GC length bins make the length-agnostic GC correction too noisy.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "0.0",
            value_parser = parse_gc_length_trim_rare,
            help_heading = "GC Correction"
        )
    )]
    pub gc_length_trim_rare: f64,

    /// 2bit reference genome file [path]
    ///
    /// NOTE: Required when specifying `--gc-file`.
    ///
    /// E.g., "hg38.2bit" from UCSC (https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit).
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

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl LengthsConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            output_prefix: String::new(),
            decimals: DEFAULT_OUTPUT_DECIMALS,
            indel_mode: IndelMode::Ignore,
            clip_mode: ClipMode::Aligned,
            max_soft_clips: DEFAULT_MAX_SOFT_CLIPS,
            max_deletion_bases: DEFAULT_MAX_DELETION_BASES,
            windows: DistributionWindowsArgs::default(),
            window_assignment: AssignToWindowArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            length_bins: vec!["30:1001:1".to_string()],
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            tile_size: 20000000,
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
            gc: ApplyGCArgFileOnly {
                gc_file: None,
                neutralize_invalid_gc: false,
            },
            gc_length_weighting: MarginalizeLengthsWeightingScheme::Equal,
            gc_length_range: GCLengthRange::Requested,
            gc_length_trim_rare: 0.0,
            ref_2bit: None,
            logging: LoggingArgs::default(),
        }
    }

    pub fn set_indel_mode(&mut self, mode: IndelMode) {
        self.indel_mode = mode;
    }

    pub fn set_windows(&mut self, windows: DistributionWindowsArgs) {
        self.windows = windows;
    }

    pub fn set_window_assignment(&mut self, assign: AssignToWindowArgs) {
        self.window_assignment = assign;
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

    /// Set exact one-base output bins for an inclusive length range.
    ///
    /// This is a programmatic convenience for tests and callers that want the
    /// old per-length distribution over `min_length..=max_length`.
    pub fn set_per_bp_length_bins(&mut self, min_length: u32, max_length: u32) {
        assert!(min_length <= max_length, "min length must be <= max length");
        let exclusive_end = max_length
            .checked_add(1)
            .expect("max length too large to build length-bin range");
        self.set_length_bins_spec(format!("{min_length}:{exclusive_end}:1"));
    }

    pub fn resolve_length_bins(&self) -> Result<Vec<u32>> {
        resolve_length_bin_edges(
            &self.length_bins,
            MIN_ACGT_BASES_FOR_GC_FRACTION,
            MAX_SUPPORTED_FRAGMENT_LENGTH,
        )
    }

    pub fn set_unpaired(&mut self, unpaired: UnpairedArgs) {
        self.unpaired = unpaired;
    }

    pub fn set_scaling_factors(&mut self, scaling_factors: Option<PathBuf>) {
        self.scale_genome.scaling_factors = scaling_factors;
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    pub fn set_decimals(&mut self, decimals: u8) {
        self.decimals = decimals;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgFileOnly) {
        self.gc = gc;
    }

    pub fn set_gc_length_weighting(
        &mut self,
        gc_length_weighting: MarginalizeLengthsWeightingScheme,
    ) {
        self.gc_length_weighting = gc_length_weighting;
    }

    pub fn set_gc_length_range(&mut self, gc_length_range: GCLengthRange) {
        self.gc_length_range = gc_length_range;
    }

    pub fn set_gc_length_trim_rare(&mut self, gc_length_trim_rare: f64) {
        self.gc_length_trim_rare = gc_length_trim_rare;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}

fn parse_gc_length_trim_rare(raw_value: &str) -> std::result::Result<f64, String> {
    let value = raw_value
        .parse::<f64>()
        .map_err(|_| "--gc-length-trim-rare must be a number".to_string())?;
    validate_gc_length_trim_rare(value).map_err(|error| error.to_string())?;
    Ok(value)
}

pub(crate) fn validate_gc_length_trim_rare(value: f64) -> Result<()> {
    anyhow::ensure!(
        value.is_finite() && (0.0..1.0).contains(&value),
        "--gc-length-trim-rare must be finite and within [0, 1)"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("config_tests.rs");
}
