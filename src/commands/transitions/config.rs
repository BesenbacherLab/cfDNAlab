use std::path::PathBuf;

use crate::{
    ToCliCommand,
    cli_command::helpers::*,
    commands::{
        cli_common::{
            ChromosomeArgs, FragmentLengthArgs, FragmentPositionSelectionArgs, IOCArgs,
            Ref2BitRequiredArgs, ScaleGenomeArgs, WindowsArgs,
        },
        fragment_kmers::config::{FragmentKmersSharedArgs, push_fragment_kmers_shared_cli_args},
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
/// Use `cfdna visualize-positions` to check what bases are counted at with various position selection settings.
///
/// ## Example
///
/// ```rust,ignore
/// // First-order transition probabilities in the 10 first bases from each 5'
/// // NOTE: To reproduce the features in `Ji et al. 2025` (https://doi.org/10.1101/2025.09.09.25335450),
/// // you can calculate the **initial probabilities** from the first transition dimension
/// // by summing all motifs that share the same set of K-1 motifs (e.g. freq_AA = sum(AAA, AAT, AAC, AAG)).
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
#[derive(Debug, Clone, PartialEq)]
pub struct TransitionsConfig {
    /// Args shared with `fragment-kmers`
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub shared_args: FragmentKmersSharedArgs,

    /// List of transition orders `[integer]`
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
}

impl TransitionsConfig {
    pub fn new(ioc: IOCArgs, ref_genome: Ref2BitRequiredArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            shared_args: FragmentKmersSharedArgs::new(
                ioc,
                ref_genome,
                chromosomes,
                "transitions".to_string(),
            ),
            orders: vec![1u8, 2u8],
            canonical: false,
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

    pub fn set_orders(&mut self, orders: Vec<u8>) {
        self.orders = orders;
    }

    pub fn set_position_selection(&mut self, position_selection: FragmentPositionSelectionArgs) {
        self.shared_args.set_position_selection(position_selection);
    }

    pub fn set_ignore_gap(&mut self, ignore_gap: bool) {
        self.shared_args.set_ignore_gap(ignore_gap);
    }

    pub fn set_canonical(&mut self, canonical: bool) {
        self.canonical = canonical;
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
}

impl ToCliCommand for TransitionsConfig {
    fn to_cli_args(&self) -> crate::Result<Vec<std::ffi::OsString>> {
        let mut args = command_args("transitions");
        push_fragment_kmers_shared_cli_args(&mut args, &self.shared_args);
        push_values(&mut args, "--orders", &self.orders);
        push_bool(&mut args, "--canonical", self.canonical);
        push_bool(&mut args, "--save-sparse", self.save_sparse);
        Ok(args)
    }
}
