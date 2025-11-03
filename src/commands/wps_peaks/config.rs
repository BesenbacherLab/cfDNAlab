use crate::commands::cli_common::ScaleGenomeArgs;
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs, WindowsArgs};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::commands::wps::config::WPSSharedConfig;

/*
What do we actually want?

WPS: Scores per position etc. Similar to coverage but in a way that allows calculating metrics yourself, etc. Perhaps bigwig files are better for that? Check for libs and size of such a file.

Peaks: Positions and stats? Just always give everything? Well, unique-positions, index-positions, stats, then allow setting multiple? For by-size, probably stats is the one to use?

*/

/// Detect nucleosome peaks via windowed protection scores (WPS) across the genome.
///
/// Only paired-end fragments with both reads present are considered.
///
/// NOTE: To extract just the WPS, see `cfdna wps` instead.
///
/// WPS: Number of fragments fully overlapping the window, minus the number of fragments ending strictly inside the window.
/// Fragments that both start and end at the exact window edges are considered fully overlapping.
///
/// ## Windowing (by-bed or by-size)
///
/// When specifying genomic windows via `--by-bed` or `--by-size`, one of the following outputs
/// is possible:
///
///  - Get the positional WPS for the included windows only (`--by-bed` *only*).
///    Excludes all positions that do not overlap a window from the output.
///    Choose between:
///     1) Indexed: Adds the original window index as an output column and keeps duplicate positions.
///     2) Unique: Overlapping windows are merged to avoid duplicate positions.
///
/// - Get the average or total WPS per window.
///
/// Without windowing, positional WPS are outputted for the selected chromosomes.
///
/// ## Blacklisting
///
/// Positions where the `--window_size` window overlaps a (dilated) blacklisted region are set to `f32::NaN` (and thus not included in sums or averages).
///
/// **Dilation**: We want to avoid any WPS scores being biased by neighbouring blacklisted intervals,
/// which can have an unreasonably high number of overlapping fragments.
/// Hence, we increase all blacklist intervals by the maximum fragment length + half the `--window_size` on both sides.
///
/// ## Scaling
///
/// When `--scaling-factors` are provided, we scale the **final per-position WPS** by the factor assigned
/// to the centre base of that position.
///
/// ## Temporary files
///
/// We write temporary files to a `<output-dir>/tmp.<output-prefix>.<random>` directory to reduce memory.
/// This directory is deleted at the end of the run. If the software is disrupted, the directory
/// may be left behind.
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
///
/// ## Examples
///
/// ```rust,ignore
///
/// // Extract peaks (these arguments are always specified, hence `...` below)
/// cfdna wps-peaks --bam <> --output-dir <> -n-threads <>
///
/// // Extract peaks in windows
/// cfdna wps-peaks ... --by-bed <> --per-window "unique-positions"
///
/// // Extract statistics about the peaks (e.g., inter-peak distances) in windows
/// cfdna wps-peaks ... --by-bed <> --per-window "stats"
///
/// ```
///
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct WPSPeaksConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub shared_args: WPSSharedConfig,

    // TODO: Make anti-dependent on skip-peaks?
    // TODO: Allow setting multiple of these? E.g. both count and average-distance? Or just "stats" for multiple peak stats? Look up what metrics people often use?
    /// What to return for peaks per window `[string]`
    ///
    /// Possible values:
    ///
    /// - `"unique-positions"`: Get the peak positions for the included windows only (`--by-bed` *only*).
    ///   Overlapping windows are merged to avoid duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"indexed-positions"`: Get the peak positions for the included windows only (`--by-bed` *only*).
    ///   Adds the original window index as an output column and keeps duplicate positions.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"count"`: Get the number of peaks per window.
    ///
    /// - `"average-distance"`: Get the average distance between peaks per window. Windows with less than 2 peaks will get a `f32::NaN` distance.
    ///
    /// **NOTE**: Ignored when no windows are specified.
    /// Required when either `--by-bed` or `--by-size` are provided and `--skip-peaks` is NOT specified.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, ignore_case = true, help_heading = "Core")
    )]
    pub per_window: Option<CoverageWindowAction>,
}

impl WPSPeaksConfig {
    pub fn new(
        ioc: IOCArgs,
        chromosomes: ChromosomeArgs,
        per_window: Option<CoverageWindowAction>,
    ) -> Self {
        Self {
            shared_args: WPSSharedConfig::new(ioc, chromosomes, "wps_peaks"),
            per_window: per_window,
        }
    }

    pub fn set_window_size(&mut self, window_size: u32) {
        self.shared_args.set_window_size(window_size);
    }

    pub fn set_decimals(&mut self, decimals: u8) {
        self.shared_args.set_decimals(decimals);
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.shared_args.set_tile_size(tile_size);
    }

    pub fn set_per_window(&mut self, action: Option<CoverageWindowAction>) {
        self.per_window = action;
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
        self.shared_args.set_windows(windows);
    }

    pub fn set_min_fragment_length(&mut self, min_fragment_length: u32) {
        self.shared_args
            .set_min_fragment_length(min_fragment_length);
    }

    pub fn set_max_fragment_length(&mut self, max_fragment_length: u32) {
        self.shared_args
            .set_max_fragment_length(max_fragment_length);
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.shared_args.set_min_mapq(min_mapq);
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.shared_args.set_require_proper_pair(require);
    }

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.shared_args.set_scale_genome(scale);
    }
}
