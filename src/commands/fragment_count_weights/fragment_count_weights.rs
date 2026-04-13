use crate::commands::{
    coverage_weights::coverage_weights::run_with_fcoverage,
    fragment_count_weights::config::FragmentCountWeightsConfig,
};
use anyhow::Result;

pub fn run(opt: &FragmentCountWeightsConfig) -> Result<()> {
    run_with_fcoverage(&opt.shared, true, "fragment_counts.scaling_factors.tsv")
}
