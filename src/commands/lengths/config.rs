use crate::{
    commands::{
        cli_common::{
            ApplyGCArgFileOnly, AssignToWindowArgs, ChromosomeArgs, FragmentLengthArgs, IOCArgs,
            ScaleGenomeArgs, UnpairedArgs, WindowsArgs,
        },
        gc_bias::correct::MarginalizeLengthsWeightingScheme,
    },
    shared::{blacklist::BlacklistStrategy, indel_mode::IndelMode},
};
use std::path::PathBuf;

/// Count fragment lengths in a BAM-file.
///
/// Fragment length: For **paired-end** sequencing, the length is defined as `end(reverse) - start(forward)`.
/// For **unpaired** sequencing where each read is a fragment, the length is defined as `[read.pos, read.end)`.
///
/// The default for windows is to count fragments by their overlap fraction. That is, most
/// fragments are counted as `1.0`, while fragments overlapping the edge of a window are counted
/// as the fraction it overlaps the window (`< 1.0`). For consequtive non-overlapping windows,
/// this conserves the total mass, as an edge-overlapping fragment will count `f` in one window
/// and `1-f` in the other window. To get base-weighted counts (i.e. coverage in the window),
/// you can multiply the output counts by their lengths (`C'[L] = L * C[L]`). **Other options**
/// include counting the full fragment if the *fragment midpoint* or a given *proportion* of
/// positions overlaps the window.
///
/// ## GC correction
///
/// Weight the contribution of each fragment based on their GC contents.
///
/// Note: The GC percentage is calculated from the full genomic coordinates (does not consider `--indel_mode`).
///
/// The length-dimension of the original correction matrix is averaged out with
/// a specifiable weighting scheme (`--gc-length-weighting`).
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

    /// How to handle insertions and deletions in fragments `[string]`
    ///
    /// Deletions: Both 'D' and 'N' in the cigar string are considered deletions.
    ///
    /// Possible values:
    ///
    /// - `"ignore"`:
    ///   Ignore whether indels are present or not. Lengths are calculated from the reference coordinates `end(reverse) - start(forward)`.
    ///
    /// - `"adjust"`:
    ///   Adjust the reference length by the observed insertions and deletions in the
    ///   observed bases (we cannot adjust in the mate-gap).
    ///   Outside the mate-overlap, all indels and deletions are adjusted for.
    ///   **Overlap**: In the mate-overlap, both reads must agree on the position-level,
    ///   with the shortest insertion selected per position.
    ///   Only overlap-positions were both reads have the indel are adjusted for.
    ///   **NOTE**: Blacklist exclusion and calculation of scaling weights (--scaling-factors)
    ///   use the full reference span.
    ///
    /// - `"skip"`:
    ///   Skip fragments with any insertion or deletion present.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "ignore",
            ignore_case = true,
            help_heading = "Core"
        )
    )]
    pub indel_mode: IndelMode,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: WindowsArgs,

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
    /// This is NOT recommended by default as it trims the tails of the length distribution.
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

    // TODO: Pretty sure about this claim, but I'm unsure whether the binning could affect this? Check
    /// How to weight the fragment length bins when estimating the global GC bias correction `[string]`
    ///
    /// To GC correct a fragment length distribution, the correction weights should be **length-agnostic**.
    /// *If* we were to use the default `fragment length bin x GC percentage bin` correction matrix, which
    /// has an independent GC correction curve per fragment length bin, the corrected counts per fragment length
    /// would have the same distributional shape as the original counts
    /// (assuming we're correcting the exact same fragments seen by `cfdna gc-bias`).
    /// We thus average out the fragment length dimension to get the overall GC bias.
    ///
    /// We have three weighting options when averaging the fragment-length-wise correction curves:
    ///     
    /// - `"equal"` weighting (default): Each fragment length bin counts the same.
    ///   This keeps the correction independent of the count distribution we're trying to estimate,
    ///   but very rare fragment length bins contribute the same as the most present fragment lengths.
    ///   For low-coverage BAM files, this *could* make the correction more volatile to outliers.
    ///
    /// - `"coverage"`-based weighting: Each fragment length bin is weighted by how often it was observed
    ///   in `cfdna gc-bias`. This should work better for the majority of the observed fragments
    ///   **but biases** the correction based on the fragment length distribution we are trying to estimate
    ///   (assuming the same BAM-file was used to estimate the GC bias).
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

impl LengthsConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            indel_mode: IndelMode::Ignore,
            windows: WindowsArgs::default(),
            window_assignment: AssignToWindowArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs::default(),
            unpaired: UnpairedArgs { reads_are_fragments: false },
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
            gc: ApplyGCArgFileOnly {
                gc_file: None,
                drop_invalid_gc: false,
            },
            gc_length_weighting: MarginalizeLengthsWeightingScheme::Equal,
            ref_2bit: None,
        }
    }

    pub fn set_indel_mode(&mut self, mode: IndelMode) {
        self.indel_mode = mode;
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
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
