use std::ops::{Deref, DerefMut};

use crate::commands::coverage_weights::scaling_weights_config::ScalingWeightsArgs;

/// Extract fragment coverage-based smoothing weights in large genomic bins ("megabins")
/// with a rolling window and calculate normalizing scaling factors for smoothing
/// the genome.
///
/// Use this when you want the smoothing profile to reflect total fragment coverage
/// across the genome. This is the natural choice for coverage-like analyses.
/// In contrast, fragment count-based weights (see `cfdna fragment-count-weights`)
/// make long and short fragments contribute more equally.
///
/// Outputs scaling factors per stride to allow other methods to apply the normalization.
/// Files are written as:
///
/// `<prefix>.coverage.scaling_factors.tsv`
///
/// **Multipliers**: After normalization of the non-zero smoothed coverage values to
/// a global mean of `1.0`, the values are **inverted** to **multiplicative** scaling factors.
///
/// ## Coverage
///
/// Internally, this command runs `fcoverage --by-size <stride> --per-window average`
/// and then smooths those stride-bin averages.
///
/// Fragment counting therefore follows `fcoverage`.
/// By default, the full fragment span is counted, except for deletions and skipped
/// regions that are not covered by the other read.
///
/// ## Fragment span definition
///
/// **Paired-end**: `[forward.pos, reverse.end)`.
///
/// **Unpaired** where each read is a fragment: `[read.pos, read.end)`.
///
/// ## GC correction
///
/// When downstream tools should use both genomic smoothing and GC-bias correction,
/// you can build the smoothing weight off GC-corrected fragment coverage by supplying either
/// `--gc-file` or `--gc-tag`. This avoids over-correction where the genomic smoothing scalars
/// partly reflect large-scale GC bias.
///
/// The written TSV records whether GC correction was used so downstream commands can check
/// whether the two transformations are used together consistently or not.
///
/// ## Smoothing
///
/// Smoothing is performed as a triangular moving average, calculating
/// a weighted average of coverages from all bins overlapping a stride.
///
/// ### Example
///
/// Assuming a bin-size of 6 and stride size of 2 (normally defaults to 5Mb and 0.5Mb respectively).
///
/// **Stride bins** (fixed along genome, each with an average coverage):
///
/// `[A] [B] [C] [D] [E] [F] [G] ...`
///
/// **Overlapping megabins** (`MB*`) (each covers 3 stride-bins). **`W_D`**, the number of overlapping megabins,
/// is the (unnormalized) weight of each stride-bin in the weighted-average coverage for stride-bin `D`:
///
/// ```text
///
/// MB1: [A][B][C]
///
/// MB2:    [B][C][D]
///
/// MB3:       [C][D][E]
///
/// MB4:          [D][E][F]
///
/// MB5:             [E][F][G]
///
/// W_D: [0][1][2][3][2][1][0]
///
/// ```
///
/// At chromosome edges, the weights are truncated (e.g., `W_D: [2][3][2][1][0]`).
///
/// The weights are normalized by their sum (after potential truncation at edges).
///
/// Stride bins with undefined average coverage, for example fully blacklisted bins from
/// `fcoverage`, are skipped while smoothing neighboring bins. They may still get a finite
/// smoothed value from neighboring support, but their scaling factor is written as `0`.
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
pub struct CoverageWeightsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub shared: ScalingWeightsArgs,

    /// Ignore inter-mate gap when building coverage-based smoothing weights `[flag]`
    ///
    /// Use this when downstream coverage analyses will also use `fcoverage --ignore-gap`.
    ///
    /// Cannot be used with `--reads-are-fragments`.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub ignore_gap: bool,
}

impl CoverageWeightsConfig {
    pub fn new(
        ioc: crate::commands::cli_common::IOCArgs,
        chromosomes: crate::commands::cli_common::ChromosomeArgs,
    ) -> Self {
        Self {
            shared: ScalingWeightsArgs::new(ioc, chromosomes),
            ignore_gap: false,
        }
    }

    pub fn set_ignore_gap(&mut self, ignore_gap: bool) {
        self.ignore_gap = ignore_gap;
    }
}

impl Deref for CoverageWeightsConfig {
    type Target = ScalingWeightsArgs;

    fn deref(&self) -> &Self::Target {
        &self.shared
    }
}

impl DerefMut for CoverageWeightsConfig {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.shared
    }
}
