use crate::{
    commands::{
        cli_common::{
            ApplyGCArgs, ChromosomeArgs, FragmentLengthArgs, IOCArgs, ScaleGenomeArgs,
            UnpairedArgs, WindowsArgs,
        },
        ends::config_structs::*,
    },
    shared::{blacklist::BlacklistStrategy, indel_mode::IndelMotifFilterPolicy},
};
use std::path::PathBuf;

/// Count fragment end motifs in a BAM-file.
///
/// Writes either:
///
/// - a dense `.npy` matrix with shape `(# windows, # motifs)` when `--all-motifs` is enabled
/// - or a sparse `.npz` matrix otherwise
///
/// along with a text file with the matching motif labels.
///
/// ## GC correction
///
/// Weight the contribution of each fragment based on their GC contents per fragment length.
///
/// ## Genomic smoothing (--scaling-factors)
///
/// Weight how genomic regions contribute to the count distribution(s), e.g., to reduce the
/// influence of copy number alterations (if that is meaningful to your analysis).
/// This weights the contribution of each fragment by region-wise precomputed scaling factors.
///
/// Can be precomputed with `cfdna coverage-weights`.
///
/// ## Window assignment
///
/// By default, a motif is counted in the window the fragment end falls in with the weight 1.0 (before correction/scaling).
///
/// Alternatively, we can weight the motif by how much the fragment overlaps the window or
/// we can count both end motifs of a fragment if the *fragment midpoint* or a given
/// *proportion* of positions overlaps the window.
///
/// ## Blacklisting
///
/// Ignores fragments that overlap blacklisted regions with a given proportion.
///
/// Motifs overlapping blacklisted regions are skipped.
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
pub struct EndsConfig {
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
    ///   `<prefix>.end_motifs.npy`
    ///   `<prefix>.end_motifs.sparse.npz`
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// 2bit reference genome file [path]
    ///
    /// NOTE: Required when using reference bases or specifying `--gc-file`.
    ///
    /// E.g., "hg38.2bit" from UCSC (https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit).
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'r',
            long,
            value_parser,
            required = false,
            help_heading = "Core"
        )
    )]
    pub ref_2bit: Option<PathBuf>,

    /// Number of bases to use from within the fragment `[integer]`
    #[cfg_attr(feature = "cli", clap(long, required = true, help_heading = "Motifs"))]
    pub k_within: usize,

    /// Number of bases to use from outside the fragment `[integer]`
    #[cfg_attr(feature = "cli", clap(long, required = true, help_heading = "Motifs"))]
    pub k_outside: usize,

    /// Whether to get the within-fragment bases from the read or the reference `[string]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "read",
            requires_if("reference", "ref_2bit"),
            help_heading = "Motifs"
        )
    )]
    pub source_within: KmerSource,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub clip: ClippingArgs,

    /// When to filter motifs due to indels.
    ///
    /// Deletions: Both 'D' and 'N' in the cigar string are considered deletions.
    ///
    /// Possible values:
    ///
    /// - `"auto"`:
    ///   Select the option based on the source.
    ///   
    ///   For **read**-sequence bases, allow indels in the alignment.
    ///
    ///   For **reference** bases, skip motifs with indels in the alignment.
    ///
    /// - `"skip-affected-end"`:
    /// Always skip motifs overlapping indels.
    ///
    /// - `"skip-affected-fragment"`:
    /// Skip **fragments** when either end overlap indels.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "auto",
            ignore_case = true,
            help_heading = "Filtering"
        )
    )]
    pub indel_filter: IndelMotifFilterPolicy,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub window_assignment: AssignMotifToWindowArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    /// Collapse each motif with its reverse-complement [flag]
    ///
    /// How:
    ///
    /// - (Always) Motifs are oriented so they run from the fragment end inward in 5'->3' direction.
    ///
    /// - Motifs are collapsed with their complement, using the lexicographically smaller motif representation.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Motifs"))]
    pub collapse_complement: bool,

    /// Include every possible motif in the output, even if its count is zero  [flag]
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Motifs"))]
    pub all_motifs: bool,

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
    /// **NOTE**: Motifs overlapping blacklisted regions are always skipped. This strategy is
    /// for further filtering of the full fragments. This is useful when you generally don't trust
    /// the reference sequences in blacklisted regions.
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
}

impl EndsConfig {
    pub fn new(
        ioc: IOCArgs,
        chromosomes: ChromosomeArgs,
        k_within: usize,
        k_outside: usize,
    ) -> Self {
        Self {
            ioc,
            output_prefix: String::new(),
            ref_2bit: None,
            k_within,
            k_outside,
            source_within: KmerSource::Read,
            clip: ClippingArgs {
                clip_strategy: ClipStrategy::Raw,
                max_soft_clips: None,
            },
            indel_filter: IndelMotifFilterPolicy::Auto,
            all_motifs: false,
            collapse_complement: false,
            windows: WindowsArgs::default(),
            window_assignment: AssignMotifToWindowArgs::default(),
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
            gc: ApplyGCArgs {
                gc_file: None,
                gc_tag: None,
                drop_invalid_gc: false,
            },
        }
    }

    pub fn set_indel_filter(&mut self, filter: IndelMotifFilterPolicy) {
        self.indel_filter = filter;
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
        self.windows = windows;
    }

    pub fn set_window_assignment(&mut self, assign: AssignMotifToWindowArgs) {
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

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}
