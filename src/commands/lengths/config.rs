use crate::{
    commands::cli_common::{
        AssignToWindowArgs, ChromosomeArgs, FragmentLengthArgs, IOCArgs, ScaleGenomeArgs,
        WindowsArgs,
    },
    shared::{blacklist::BlacklistStrategy, indel_mode::IndelMode},
};
use std::path::PathBuf;

/// Count fragment lengths in a BAM-file.
///
/// Fragment length is defined as `end(reverse) - start(forward)`.
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
pub struct LengthsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

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
            default_value = "reference",
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

impl LengthsConfig {
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            indel_mode: IndelMode::Ignore,
            windows: WindowsArgs::default(),
            window_assignment: AssignToWindowArgs::default(),
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

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }
}
