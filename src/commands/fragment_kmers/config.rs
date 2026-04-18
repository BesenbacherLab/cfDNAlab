use std::path::PathBuf;

use crate::{
    commands::cli_common::{
        ApplyGCArgs, BaseSelectionArgs, ChromosomeArgs, FragmentLengthArgs,
        FragmentPositionSelectionArgs, IOCArgs, LoggingArgs, Ref2BitRequiredArgs, ScaleGenomeArgs,
        UnpairedArgs, WindowsArgs,
    },
    shared::{
        blacklist::BlacklistStrategy,
        indel_mode::IndelMode,
        positioning::{BasesFrom, MismatchBasesFrom, ReferenceFrame},
    },
};

/// Commands that are shared with `transitions``
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct FragmentKmersSharedArgs {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ref_genome: Ref2BitRequiredArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.k3_counts.npy`,
    ///   `<prefix>.k3_motifs.txt`,
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'x',
            default_value_t = String::new(),
            hide_default_value = true,
            help_heading = "Core"
        )
    )]
    pub output_prefix: String,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    // TODO: Align with `cfdna ends` ("adjust" makes no sense here)
    // TODO: Is it still correct that scaling weights use the full reference span? 5th last line
    /// How to handle insertions and deletions in fragments `[string]`
    ///
    /// Deletions: Both 'D' and 'N' in the cigar string are considered deletions.
    ///
    /// Possible values:
    ///
    /// - `"ignore"`:
    ///   Ignore whether indels are present or not.
    ///   Kmers are extracted for the full/offset fragment span from the reference genome.
    ///
    /// - `"adjust"`:
    ///   Adjust the counts by excluding k-mers overlapping positions with observed insertions and deletions in the
    ///   observed bases (we cannot adjust in mate-gaps).
    ///   Outside the mate-overlap, all indels and deletions are adjusted for.
    ///   **Overlap**: In the mate-overlap, both reads must agree on the position-level.
    ///   Only overlap-positions were both reads have the indel are excluded.
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

    /// Ignore inter-mate gap `[flag]`
    ///
    /// Disable counting in the gap between reads (i.e., `[forward.end, reverse.start)`)
    /// when the two reads do not overlap.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub ignore_gap: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub position_selection: FragmentPositionSelectionArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub base_selection: BaseSelectionArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    // TODO: Add that we use the scaling weight for the first kmer-position
    // And that sf=0 for any kmer base guarantees the kmer is excluded
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
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// Two levels of filtering are performed. First, all blacklisted regions are assigned
    /// the N-"base" to exclude k-mers that include the positions. Then, depending on the `--blacklist-strategy`,
    /// fragments overlapping blacklisted regions with some fraction are excluded.
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

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl FragmentKmersSharedArgs {
    pub fn new(
        ioc: IOCArgs,
        ref_genome: Ref2BitRequiredArgs,
        chromosomes: ChromosomeArgs,
        output_prefix: String,
    ) -> Self {
        Self {
            ioc,
            ref_genome,
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            output_prefix,
            tile_size: 20_000_000,
            position_selection: FragmentPositionSelectionArgs {
                frame: vec![ReferenceFrame::Left],
                positions: vec!["..".to_string()],
                step: vec![1],
            },
            base_selection: BaseSelectionArgs {
                bases_from: BasesFrom::Reference,
                mismatch_bases_from: MismatchBasesFrom::NearestRead,
            },
            indel_mode: IndelMode::Ignore,
            ignore_gap: false,
            windows: WindowsArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs::default(),
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
            logging: LoggingArgs::default(),
        }
    }

    pub fn set_output_prefix(&mut self, output_prefix: String) {
        self.output_prefix = output_prefix;
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

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    pub fn set_position_selection(&mut self, position_selection: FragmentPositionSelectionArgs) {
        self.position_selection = position_selection;
    }

    pub fn set_base_selection(&mut self, base_selection: BaseSelectionArgs) {
        self.base_selection = base_selection;
    }

    pub fn set_ignore_gap(&mut self, ignore_gap: bool) {
        self.ignore_gap = ignore_gap;
    }

    pub fn set_indel_mode(&mut self, indel_mode: IndelMode) {
        self.indel_mode = indel_mode;
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
        self.windows = windows;
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.scale_genome = scale;
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

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }
}

// TODO: Add minimum or min mean base quality filtering!

/// Count k-mers within the fragments in a BAM-file.
///
/// **Experimental**: enable via `--features cmd_fragment_kmers` during `cargo build/install`.
///
/// Whereas the `cfdna ends` tool extracts end-motifs, this tool extracts k-mers
/// from all the **selected positions** in the fragment.
///
/// Use `cfdna visualize-positions` to check what bases are counted at with various position selection settings.
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
pub struct FragmentKmersConfig {
    /// Args shared with downstream tools like `transitions`
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub shared_args: FragmentKmersSharedArgs,

    /// List of K-mer sizes [integer].
    ///
    /// When counting for many kmer-sizes (>8), consider splitting
    /// into multiple runs to reduce memory consumption at a time.
    ///
    /// Example: `--kmer-sizes 3 5 11`
    #[cfg_attr(
        feature = "cli",
        clap(short = 'k', long, num_args = 1.., value_parser = clap::value_parser!(u8).range(1..28), required=true, help_heading="Core"))]
    pub kmer_sizes: Vec<u8>,

    // TODO: Re-audit whether reverse-complement collapse is the right biological contract here.
    // `ends` intentionally decodes to a final orientation first and then compares against the
    // same-orientation complement, but `fragment-kmers` may still want true reverse-complement
    // collapse. Decide that explicitly before adding more contract-level tests or docs here.
    /// Collapse each kmer with its reverse-complement. [flag]
    ///
    /// Each k-mer is compared with its reverse complement and the lexicographically smaller
    /// motif is kept.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub canonical: bool,

    /// Enable positional counting based on --frame/--positions `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub positional_counts: bool,

    /// Save counts as sparse-array. [flag]
    ///
    /// For large kmer-sizes, we cannot save dense arrays with all motifs
    /// unless we have a LOT of RAM and storage space. Enable this
    /// flag to save as a COO sparse array that can be opened in
    /// python via `scipy.sparse.load_npz()`.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub save_sparse: bool,
    // #[cfg_attr(feature = "cli", clap(flatten))]
    // gc: GCArgs,

    // #[cfg_attr(feature = "cli", clap(flatten))]
    // two_bit: TwoBitArgs,
}

impl FragmentKmersConfig {
    pub fn new(ioc: IOCArgs, ref_genome: Ref2BitRequiredArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            shared_args: FragmentKmersSharedArgs::new(
                ioc,
                ref_genome,
                chromosomes,
                "fragment_kmers".to_string(),
            ),
            kmer_sizes: vec![3u8],
            canonical: false,
            positional_counts: false,
            save_sparse: false,
        }
    }

    pub fn set_output_prefix(&mut self, output_prefix: String) {
        self.shared_args.set_output_prefix(output_prefix);
    }

    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.shared_args.set_blacklist(blacklist);
    }

    pub fn set_blacklist_min_size(&mut self, blacklist_min_size: u64) {
        self.shared_args.set_blacklist_min_size(blacklist_min_size);
    }

    pub fn set_blacklist_strategy(&mut self, blacklist_strategy: BlacklistStrategy) {
        self.shared_args.set_blacklist_strategy(blacklist_strategy);
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.shared_args.set_tile_size(tile_size);
    }

    pub fn set_kmer_sizes(&mut self, kmer_sizes: Vec<u8>) {
        self.kmer_sizes = kmer_sizes;
    }

    pub fn set_position_selection(&mut self, position_selection: FragmentPositionSelectionArgs) {
        self.shared_args.set_position_selection(position_selection);
    }

    pub fn set_base_selection(&mut self, base_selection: BaseSelectionArgs) {
        self.shared_args.set_base_selection(base_selection);
    }

    pub fn set_ignore_gap(&mut self, ignore_gap: bool) {
        self.shared_args.set_ignore_gap(ignore_gap);
    }

    pub fn set_canonical(&mut self, canonical: bool) {
        self.canonical = canonical;
    }

    pub fn set_positional_counts(&mut self, positional: bool) {
        self.positional_counts = positional;
    }

    pub fn set_save_sparse(&mut self, save_sparse: bool) {
        self.save_sparse = save_sparse;
    }

    pub fn set_indel_mode(&mut self, indel_mode: IndelMode) {
        self.shared_args.set_indel_mode(indel_mode);
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
        self.shared_args.set_windows(windows);
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.shared_args.set_scale_genome(scale);
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        self.shared_args.fragment_lengths_mut()
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.shared_args.set_min_mapq(min_mapq);
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.shared_args.set_require_proper_pair(require);
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.shared_args.set_gc(gc);
    }
}
