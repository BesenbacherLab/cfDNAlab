use crate::commands::coverage_weights::{
    config::push_scaling_weights_cli_args, scaling_weights_config::ScalingWeightsArgs,
};
use crate::{ToCliCommand, cli_command::helpers::*};
use std::ops::{Deref, DerefMut};

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
/// **Multipliers**: After normalization of the non-zero smoothed fragment counts to
/// a global mean of `1.0`, the values are inverted to multiplicative scaling factors.
///
/// ## Fragment counts
///
/// Internally, this command runs:
///
/// `fcoverage --normalize-by-length=unit-mass --by-size <stride> --per-window total`
///
/// and then smooths those stride-bin totals.
///
/// The resulting stride-bin values are **fractional** fragment counts in each stride bin.
/// Each fragment contributes `1.0` in total. If it crosses a stride-bin boundary,
/// that contribution is split between the bins as the fractional overlap of each bin.
///
/// ## Fragment span definition
///
/// **Paired-end**: `[forward.pos, reverse.reference_end)`, the reference span
/// from the first aligned position on the forward read to the last aligned
/// position on the reverse read.
///
/// **Unpaired** where each read is a fragment: `[read.pos, read.reference_end)`,
/// the reference span from the first to the last aligned position on the read.
///
/// ## GC correction
///
/// When downstream tools should use both genomic smoothing and GC-bias correction,
/// supply `--gc-file` or `--gc-tag` here too. The command then uses corrected fragment
/// counts, which avoids over-correction downstream when the genomic smoothing factors
/// partly reflect large-scale GC bias.
///
/// The written TSV records whether GC correction was used so downstream commands can check
/// whether the two transformations are used together consistently or not.
///
/// ## Smoothing
///
/// Smoothing is performed as a triangular moving average, calculating
/// a weighted average of fragment counts from all bins overlapping a stride.
///
/// ### Example
///
/// Assuming a bin-size of 6 and stride size of 2 (normally defaults to 5Mb and 0.5Mb respectively).
///
/// **Stride bins** (fixed along genome, each with a fragment count):
///
/// `[A] [B] [C] [D] [E] [F] [G] ...`
///
/// **Overlapping megabins** (`MB*`) (each covers 3 stride-bins). **`W_D`**, the number of overlapping megabins,
/// is the (unnormalized) weight of each stride-bin in the weighted-average fragment count for stride-bin `D`:
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
/// The stride bins are further weighted by their number of eligible bases (non-blacklisted
/// positions). This also handles the often shorter final stride bin per chromosome.
///
/// The weights are normalized by their sum (after potential truncation at edges).
///
/// Fully blacklisted stride bins are skipped while smoothing neighboring bins. They may still get
/// a finite smoothed value from neighboring support, but their scaling factor is written as `0`.
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
#[derive(Debug, Clone, PartialEq)]
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

impl ToCliCommand for FragmentCountWeightsConfig {
    fn to_cli_args(&self) -> crate::Result<Vec<std::ffi::OsString>> {
        let mut args = command_args("fragment-count-weights");
        push_scaling_weights_cli_args(&mut args, &self.shared);
        Ok(args)
    }
}
