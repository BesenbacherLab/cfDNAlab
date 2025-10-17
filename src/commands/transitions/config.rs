use std::path::PathBuf;

use crate::{
    commands::{
        cli_common::{
            ChromosomeArgs, FragmentLengthArgs, FragmentPositionSelectionArgs, IOCArgs,
            Ref2BitRequiredArgs, ScaleGenomeArgs, WindowsArgs,
        },
        visualize_positions::{BasesFrom, MismatchBasesFrom, ReferenceFrame},
    },
    shared::{blacklist::BlacklistStrategy, indel_mode::IndelMode},
};

/// Calculate positional Nth-order transition probabilities within the fragment in a BAM-file.
///
/// This command wraps `cfdna fragment-kmers` and calculates the probabilities based its the output.
///
/// Pipeline: A) **Count** k-mers of size `order + 1` (e.g. 2-mers for first-order transitions) per specified position.
/// B) Calculate position-wise frequencies of all k-mers.
///
/// ## Example
///
/// ```rust,ignore
/// // First-order transition probabilities in the 10 first bases from each 5' end:
/// cfdna transitions --bam <> --output-dir <> --ref-2bit <> --n-threads 12 --orders 1 --frame nearest --positions '..10' --indel-mode adjust
/// ```
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
pub struct TransitionsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ref_genome: Ref2BitRequiredArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.n1_frequencies.npy`,
    ///   `<prefix>.n1_motifs.txt`,
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'x',
            default_value = "transitions",
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

    /// List of transition orders [integer].
    ///
    /// E.g. if you want to predict based on only the previous base, it's first order (e.g. [A>C]). This practically leads to 2-mer frequencies.
    /// If you want to predict based on the previous TWO bases, it's second order (e.g. [AT>C]).
    ///
    /// When counting for many orders (>8), consider splitting
    /// into multiple runs to reduce memory consumption at a time.
    ///
    /// Example: `--orders 1 2`
    #[cfg_attr(
        feature = "cli",
        clap(short = 'n', long, num_args = 1.., default_values_t = [1u8, 2u8], value_parser = clap::value_parser!(u8).range(1..27), help_heading="Core"))]
    pub orders: Vec<u8>,

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

    // TODO: Perhaps this should collapse around the first/terminal base instead?
    /// Collapse each kmer with its reverse-complement. [flag]
    ///
    /// Odd-sized k-mers are collapsed such that the middle base is `A` or `C`.
    /// Even-sized k-mers are collapsed to the lexicographically lowest motif.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub canonical: bool,

    /// Save counts as sparse-array. [flag]
    ///
    /// For large kmer-sizes, we cannot save dense arrays with all motifs
    /// unless we have a LOT of RAM and storage space. Enable this
    /// flag to save as a COO sparse array that can be opened in
    /// python via `scipy.sparse.load_npz()`.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub save_sparse: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub position_selection: FragmentPositionSelectionArgs,

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
    /// This is NOT recommended by default as it trims the tails of the length distribution.
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
    ///     "any", "all", "midpoint", or "proportion=<threshold>"
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
    // #[cfg_attr(feature = "cli", clap(flatten))]
    // gc: GCArgs,

    // #[cfg_attr(feature = "cli", clap(flatten))]
    // two_bit: TwoBitArgs,
}

impl TransitionsConfig {
    pub fn new(ioc: IOCArgs, ref_genome: Ref2BitRequiredArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            ref_genome,
            output_prefix: "fragment_kmers".to_string(),
            tile_size: 20_000_000,
            orders: vec![1u8, 2u8],
            position_selection: FragmentPositionSelectionArgs {
                frame: ReferenceFrame::Left,
                positions: "..".to_string(),
                step: 1,
                bases_from: BasesFrom::Reference,
                mismatch_bases_from: MismatchBasesFrom::NearestRead,
            },
            indel_mode: IndelMode::Ignore,
            ignore_gap: false,
            canonical: false,
            save_sparse: false,
            windows: WindowsArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs {
                min_fragment_length: 20,
                max_fragment_length: 1000,
            },
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
        }
    }

    pub fn set_output_prefix(&mut self, output_prefix: String) {
        self.output_prefix = output_prefix;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    pub fn set_orders(&mut self, orders: Vec<u8>) {
        self.orders = orders;
    }

    pub fn set_position_selection(&mut self, position_selection: FragmentPositionSelectionArgs) {
        self.position_selection = position_selection;
    }

    pub fn set_ignore_gap(&mut self, ignore_gap: bool) {
        self.ignore_gap = ignore_gap;
    }

    pub fn set_canonical(&mut self, canonical: bool) {
        self.canonical = canonical;
    }

    pub fn set_save_sparse(&mut self, save_sparse: bool) {
        self.save_sparse = save_sparse;
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
}
