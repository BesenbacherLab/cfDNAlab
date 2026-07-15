#[test]
#[ignore = "requires CFDNALAB_BENCH_FCOVERAGE_TSV and measures wall-clock performance"]
fn benchmark_background_reading_fcoverage_output_loader() -> Result<()> {
    let path =
        crate::background_reading_benchmarks::required_path("CFDNALAB_BENCH_FCOVERAGE_TSV")?;
    let group_index_path = crate::background_reading_benchmarks::optional_path(
        "CFDNALAB_BENCH_FCOVERAGE_GROUP_INDEX",
    )?;
    crate::background_reading_benchmarks::compare_read_modes(
        "fcoverage output loader",
        |read_in_background| {
            FCoverageParser::new(&path, group_index_path.as_deref(), read_in_background).load()
        },
    )
}
