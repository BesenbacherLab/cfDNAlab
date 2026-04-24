use std::ops::{Deref, DerefMut};

use crate::commands::coverage_weights::scaling_weights_config::ScalingWeightsArgs;

/// Extract fragment count-based smoothing weights in large genomic bins ("megabins")
/// with a rolling window and calculate normalizing scaling factors for smoothing
/// the genome.
///
/// Use this when you want long and short fragments to contribute more equally to the
/// large-scale smoothing profile. In contrast, regular coverage-based weights
/// (see `cfdna coverage-weights`) count long fragments higher simply because
/// they cover more bases.
///
/// Outputs scaling factors per stride to allow other methods to apply the normalization.
/// Files are written as:
///
/// `<prefix>.fragment_counts.scaling_factors.tsv`
///
/// The scaling factors are *inverted*, so normalization becomes multiplication.
/// Zero-valued smoothed fragment mass leads to zero-valued scaling factors.
/// Non-zero factors have `mean == 1.0`.
///
/// ## Fragment counts
///
/// Internally, this command runs:
///
/// `fcoverage --normalize-by-length=unit-mass --by-size <stride> --per-window total`
///
/// (The `unit-mass` mode is used as it's cheaper than rescaling and normalizes to the same weights.)
///
/// and then smooths those stride-bin totals.
///
/// The resulting stride-bin values approximate fragment counts in each stride bin.
/// A full fragment contributes total mass 1.0, split across the stride bins it overlaps
/// according to covered span.
///
/// Strictly speaking this is still an approximation since fragments overlapping
/// multiple stride bins are counted partly in each, but in sufficiently large
/// bins the approximation error is tiny.
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
/// you can build the smoothing weight off GC-corrected fragment mass by supplying either
/// `--gc-file` or `--gc-tag`. This avoids over-correction where the genomic smoothing scalars
/// partly reflect large-scale GC bias.
///
/// The written TSV records whether GC correction was used so downstream commands can check
/// whether the two transformations are used together consistently or not.
///
/// ## Smoothing
///
/// Smoothing is performed as a triangular moving average, calculating
/// a weighted average of fragment-mass values from all bins overlapping a stride.
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
pub struct FragmentCountWeightsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub shared: ScalingWeightsArgs,
}

impl FragmentCountWeightsConfig {
    pub fn new(
        ioc: crate::commands::cli_common::IOCArgs,
        chromosomes: crate::commands::cli_common::ChromosomeArgs,
    ) -> Self {
        Self {
            shared: ScalingWeightsArgs::new(ioc, chromosomes),
        }
    }
}

impl Deref for FragmentCountWeightsConfig {
    type Target = ScalingWeightsArgs;

    fn deref(&self) -> &Self::Target {
        &self.shared
    }
}

impl DerefMut for FragmentCountWeightsConfig {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.shared
    }
}
