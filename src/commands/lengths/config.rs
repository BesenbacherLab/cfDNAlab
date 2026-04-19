use crate::{
    commands::{
        cli_common::{
            ApplyGCArgFileOnly, AssignToWindowArgs, ChromosomeArgs, DistributionWindowsArgs,
            FragmentLengthArgs, IOCArgs, LoggingArgs, ScaleGenomeArgs, UnpairedArgs,
        },
        gc_bias::correct::MarginalizeLengthsWeightingScheme,
    },
    shared::{blacklist::BlacklistStrategy, clip_mode::ClipMode, indel_mode::IndelMode},
};
use std::path::PathBuf;

pub const DEFAULT_MAX_SOFT_CLIPS: u16 = 256;
pub const MAX_MAX_SOFT_CLIPS: u16 = 256;

// TODO: Add length bins to enable e.g. short/long in much smaller windows without exploding RAM

/// Count fragment lengths in a BAM-file.
///
/// Writes an `.npy` file with shape (# windows, # lengths).
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
/// Note: The GC percentage is calculated from the **aligned** reference span.
/// It does not consider `--indel-mode` or `--clip-mode`.
///
/// The length-dimension of the original correction matrix is averaged out with
/// a specifiable weighting scheme (`--gc-length-weighting`).
///
/// ## Genomic smoothing (--scaling-factors)
///
/// Weight how genomic regions contribute to the length distribution(s), e.g., to reduce the
/// influence of copy number alterations. This weights the contribution of each fragment
/// by region-wise precomputed scaling factors.
///
/// Can be precomputed with `cfdna coverage-weights` or `cfdna fragment-count-weights`.
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
/// To get base-weighted counts (i.e. coverage in the window), you can multiply the output
/// counts by their lengths (`C'[L] = L * C[L]`; Remember to account for the minimum fragment
/// length offset).
///
/// Other options include counting the full fragment if the *fragment midpoint* or a given
/// *proportion* of positions overlaps the window.
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
    ///   `<prefix>.length_counts.npy`
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, help_heading = "Core")
    )]
    pub output_prefix: String,

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

    /// Skip fragments where one or both ends have more soft-clipped bases than this `[integer]`
    ///
    /// Use `--clip-mode skip` to discard all soft-clipped fragments.
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

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

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

    /// How to weight the fragment length bins when estimating the global GC bias correction `[string]`
    ///
    /// To GC correct a fragment length distribution, the correction weights should be **length-agnostic**.
    ///
    /// The default `fragment-length x GC` matrix has one correction curve per length bin,
    /// so using it would preserve the original length distribution (assuming we're correcting the
    /// same fragments seen by `cfdna gc-bias`).
    ///
    /// We therefore average out the fragment length dimension to get a single, length-agnostic GC bias curve.
    ///
    /// We have three weighting options when averaging the fragment-length-wise correction curves:
    ///     
    /// - `"equal"` weighting (default): Weight every length bin the same.
    ///
    ///   Keeps the correction independent of the distribution we are trying to estimate.
    ///
    ///   Downside: Rare fragment length bins contribute the same as the most present fragment lengths.
    ///   
    ///   For low-coverage BAM files, this *could* make the correction more volatile to outliers.
    ///
    /// - `"coverage"`-based weighting: Weight lengths by how often they were observed in `cfdna gc-bias`.
    ///
    ///   This should work better for the majority of the observed fragments **BUT**:
    ///
    ///   Downside: **Biases** the correction based on the length distribution we are trying to estimate.
    ///
    /// - `"max-coverage"` weighting: Use the GC curve for the most-observed fragment length bin.
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
            indel_mode: IndelMode::Ignore,
            clip_mode: ClipMode::Aligned,
            max_soft_clips: DEFAULT_MAX_SOFT_CLIPS,
            windows: DistributionWindowsArgs::default(),
            window_assignment: AssignToWindowArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs::default(),
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

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
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

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}
