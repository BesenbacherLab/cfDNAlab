#[test]
#[ignore = "requires CFDNALAB_BENCH_LENGTHS_TSV and measures wall-clock performance"]
fn benchmark_background_reading_lengths_output_loader() -> Result<()> {
    let path = crate::background_reading_benchmarks::required_path("CFDNALAB_BENCH_LENGTHS_TSV")?;
    crate::background_reading_benchmarks::compare_read_modes(
        "lengths output loader",
        |read_in_background| LengthsParser::new(&path, read_in_background).load(),
    )
}
