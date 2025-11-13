use crate::commands::cli_common::ScaleGenomeArgs;
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs, WindowsArgs};
use crate::commands::wps::config::WPSSharedConfig;
use crate::commands::wps_peaks::window_peak_results::PeaksWindowAction;

/*
What do we actually want?

WPS: Scores per position etc. Similar to coverage but in a way that allows calculating metrics yourself, etc. Perhaps bigwig files are better for that? Check for libs and size of such a file.

Peaks: Positions and stats? Just always give everything? Well, unique-positions, index-positions, stats, then allow setting multiple? For by-size, probably stats is the one to use?

*/

/// Detect nucleosome peaks via windowed protection scores (WPS) across the genome.
/// 
/// **Experimental**: enable via `--features cmd_wps_peaks cmd_wps` during `cargo build/install`.
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
/// ## Smoothing
///
/// The WPS values are smoothed with a Savitzky-Golar filter (second order polynomial, 21bp window filter), as used in Snyder et al.
///
/// All masked positions are edges to the smoother. [TODO: Describe this edge part (and masking, reference blacklisting)]
///
/// Disable smoothing with `--no-smoothing` if you want to keep the raw WPS values.
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

    // TODO: Allow setting multiple of these? E.g. both count and average-distance? Or just "stats" for multiple peak stats? Look up what metrics people often use?
    /// What to return for peaks per window `[string]`
    ///
    /// Possible values:
    ///
    /// - `"unique-positions"`: Get the distinct peak coordinates inside the provided windows (`--by-bed` only).
    ///   Overlapping windows are merged before processing to avoid duplicate rows.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"indexed-positions"`: Emit all peak coordinates inside the provided windows together
    ///   with the original window index (`--by-bed` only). Overlapping windows keep duplicates so each window is
    ///   reported independently.
    ///   Excludes all positions that do not overlap a window from the output.
    ///
    /// - `"stats"`: Emit peak counts as well as average and median distances between peaks per
    ///   window.
    ///
    /// Ignored when no windows are supplied. Required whenever `--by-bed` or `--by-size` are
    /// provided.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, ignore_case = true, help_heading = "Core")
    )]
    pub per_window: Option<PeaksWindowAction>,

    /// Size of window for normalizing the WPS values before smoothing `[integer]`
    ///
    /// A rolling median of this width is subtracted from the (optionally smoothed) WPS signal.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1000", value_parser = clap::value_parser!(u32).range(100..),  help_heading = "Core")
    )]
    pub normalize_bp: u32,

    /// Minimum usable bases required inside the normalization window `[integer]`
    ///
    /// Windows with fewer valid bases yield `NaN` after normalization.
    #[cfg_attr(
        feature = "cli",
        clap(long = "min-unmasked", default_value = "400", value_parser = clap::value_parser!(u32).range(1..), help_heading = "Core")
    )]
    pub min_unmasked: u32,

    /// Disable Savitzky-Golay smoothing of the WPS signal `[flag]`
    ///
    /// Smoothing is enabled by default to reproduce Snyder et al.
    #[cfg_attr(
        feature = "cli",
        clap(long = "no-smoothing", action, help_heading = "Core")
    )]
    pub no_smoothing: bool,

    // TODO: revisit default after empirical tuning
    /// Minimum residual height required to keep a peak `[float]`
    ///
    /// Any peak whose baseline-adjusted WPS height (residual) is smaller than this value
    /// is discarded. Lower values keep weaker peaks (useful for low-coverage cfDNA), higher values keep only stronger peaks.
    ///
    /// NOTE: To make this more sequencing depth agnostic, use genomic smoothing (`--scaling-factors`) to
    /// make the average coverage `~1.0` and then tune this minimum to that value. This should
    /// generalize better across samples.
    #[cfg_attr(
        feature = "cli",
        clap(
            long = "min-peak-height",
            default_value = "5.0",
            value_parser = parse_nonnegative_f32,
            help_heading = "Core"
        )
    )]
    pub min_peak_height: f32,
}

impl WPSPeaksConfig {
    pub fn new(
        ioc: IOCArgs,
        chromosomes: ChromosomeArgs,
        per_window: Option<PeaksWindowAction>,
    ) -> Self {
        Self {
            shared_args: WPSSharedConfig::new(ioc, chromosomes, "wps_peaks"),
            per_window: per_window,
            normalize_bp: 1000,
            min_unmasked: 400,
            no_smoothing: false,
            min_peak_height: 5.0,
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

    pub fn set_per_window(&mut self, action: Option<PeaksWindowAction>) {
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

    pub fn set_min_unmasked(&mut self, min_unmasked: u32) {
        self.min_unmasked = min_unmasked;
    }

    pub fn set_no_smoothing(&mut self, no_smoothing: bool) {
        self.no_smoothing = no_smoothing;
    }

    pub fn set_min_peak_height(&mut self, min_peak_height: f32) {
        self.min_peak_height = min_peak_height;
    }
}

#[cfg_attr(not(feature = "cli"), allow(dead_code))]
fn parse_nonnegative_f32(input: &str) -> Result<f32, String> {
    let value: f32 = input
        .parse()
        .map_err(|err: std::num::ParseFloatError| err.to_string())?;
    if value < 0.0 {
        Err("min-peak-height must be non-negative".into())
    } else {
        Ok(value)
    }
}
