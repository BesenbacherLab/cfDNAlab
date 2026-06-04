use crate::{
    command_run::RunOptions,
    commands::{
        coverage_weights::coverage_weights::{
            ScalingWeightsCommand, ScalingWeightsRunResult, run_with_fcoverage,
        },
        fragment_count_weights::config::FragmentCountWeightsConfig,
    },
};
use anyhow::Result;

/// Result from `fragment-count-weights`.
///
/// This is the same result shape used by `coverage-weights`, with the output representing
/// smoothed unit-fragment mass instead of average coverage.
pub type FragmentCountWeightsRunResult = ScalingWeightsRunResult;

/// Run the `fragment-count-weights` command.
///
/// This command estimates broad fragment-count structure and writes genomic scaling factors. It
/// reuses the scaling-weights implementation with length normalization enabled, so each accepted
/// fragment contributes unit mass across its span.
///
/// Reporting is controlled by `options`. `report_statistics` prints the final summary and
/// `log_statuses` controls status messages. This command does not use progress bars.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `fragment-count-weights` command.
/// - `options`:
///     Reporting controls for statistics and status logs.
///
/// Returns
/// -------
/// - `Ok(FragmentCountWeightsRunResult)`:
///     The scaling-factor path, internal `fcoverage` result, and counters.
///
/// Errors
/// ------
/// Returns an error if internal `fcoverage` fails, the intermediate TSV is malformed, or the final
/// scaling output cannot be written.
pub fn run_fragment_count_weights(
    opt: &FragmentCountWeightsConfig,
    options: RunOptions,
) -> Result<FragmentCountWeightsRunResult> {
    run_with_fcoverage(
        &opt.shared,
        true,
        ScalingWeightsCommand::FragmentCount,
        None,
        options,
    )
}
