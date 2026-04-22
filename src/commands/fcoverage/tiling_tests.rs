use super::{
    build_summary_prefixes, coverage_sum_and_counts, finalize_value, merge_scaled_positional_tiles,
};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::shared::{coverage::Coverage, interval::Interval, io::open_text_reader};
use anyhow::Result;
use std::io::{Read, Write};
use tempfile::TempDir;

fn write_zstd_lines(temp_dir: &TempDir, file_name: &str, lines: &[&str]) -> Result<()> {
    let path = temp_dir.path().join(file_name);
    let file = std::fs::File::create(path)?;
    let mut encoder = zstd::Encoder::new(file, 3)?;
    for line in lines {
        writeln!(encoder, "{line}")?;
    }
    encoder.finish()?;
    Ok(())
}

fn read_text(path: &std::path::Path) -> Result<String> {
    let mut reader = open_text_reader(path)?;
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    Ok(text)
}

fn dense_bedgraph(text: &str, chromosome: &str, chromosome_length: usize) -> Vec<f64> {
    let mut coverage = vec![0.0; chromosome_length];
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<_> = line.split('\t').collect();
        if cols[0] != chromosome {
            continue;
        }
        let start: usize = cols[1].parse().expect("start must parse");
        let end: usize = cols[2].parse().expect("end must parse");
        let value: f64 = cols[3].parse().expect("value must parse");
        for position in start..end {
            coverage[position] = value;
        }
    }
    coverage
}

#[test]
fn build_summary_prefixes_without_mask_omits_unmasked_prefixes() {
    // Arrange
    // Use a tiny finalized coverage track:
    //   index: 0   1    2   3
    //   cov:   0, 1.5, 0, 2
    //
    // Hand-derived prefixes:
    // - sum_of_squares_all:
    //     [0,
    //      0^2,
    //      0^2 + 1.5^2,
    //      0^2 + 1.5^2 + 0^2,
    //      0^2 + 1.5^2 + 0^2 + 2^2]
    //   = [0, 0, 2.25, 2.25, 6.25]
    // - nonzero_all:
    //   = [0, 0, 1, 1, 2]
    let mut coverage = Coverage::new(4);
    coverage.finalize_coverage(false);
    coverage
        .coverage_mut()
        .expect("coverage should be available after finalization")
        .copy_from_slice(&[0.0, 1.5, 0.0, 2.0]);

    // Act
    let prefixes = build_summary_prefixes(&coverage).expect("summary prefixes");

    // Assert
    assert_eq!(prefixes.sum_of_squares_all, vec![0.0, 0.0, 2.25, 2.25, 6.25]);
    assert_eq!(prefixes.nonzero_all, vec![0, 0, 1, 1, 2]);
    assert_eq!(prefixes.sum_of_squares_unmasked, None);
    assert_eq!(prefixes.nonzero_unmasked, None);
}

#[test]
fn build_summary_prefixes_with_mask_tracks_all_and_unmasked_prefixes() {
    // Arrange
    // Finalized coverage track:
    //   index: 0  1  2  3
    //   cov:   1, 0, 2, 3
    //
    // Blacklist [1, 3), so allowed positions are indices 0 and 3.
    //
    // Hand-derived prefixes:
    // - sum_of_squares_all:
    //   [0, 1, 1, 5, 14]
    // - nonzero_all:
    //   [0, 1, 1, 2, 3]
    // - sum_of_squares_unmasked:
    //   allowed values are [1, 3], so
    //   [0, 1, 1, 1, 10]
    // - nonzero_unmasked:
    //   [0, 1, 1, 1, 2]
    let mut coverage = Coverage::new(4);
    coverage.finalize_coverage(false);
    coverage
        .coverage_mut()
        .expect("coverage should be available after finalization")
        .copy_from_slice(&[1.0, 0.0, 2.0, 3.0]);
    coverage
        .set_blacklist_mask(&[Interval::new(1, 3).expect("valid blacklist interval")])
        .expect("blacklist mask");

    // Act
    let prefixes = build_summary_prefixes(&coverage).expect("summary prefixes");

    // Assert
    assert_eq!(prefixes.sum_of_squares_all, vec![0.0, 1.0, 1.0, 5.0, 14.0]);
    assert_eq!(prefixes.nonzero_all, vec![0, 1, 1, 2, 3]);
    assert_eq!(
        prefixes.sum_of_squares_unmasked,
        Some(vec![0.0, 1.0, 1.0, 1.0, 10.0])
    );
    assert_eq!(prefixes.nonzero_unmasked, Some(vec![0, 1, 1, 1, 2]));
}

#[test]
fn coverage_sum_and_counts_uses_allowed_prefixes_when_masked_indexes_exist() {
    // Arrange
    // Per-base coverage is [1, 2, 3, 4]
    // Allowed bases are positions 0 and 2, so masked coverage is [1, 0, 3, 0]
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];
    let psum_allowed = [0.0, 1.0, 1.0, 4.0, 4.0];
    let allowed_count_prefix = [0_u32, 1, 1, 2, 2];
    let blacklist_mask = [0_u8, 1, 0, 1];

    // Act
    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
        0,
        4,
        true,
        &psum_all,
        Some(&psum_allowed),
        Some(&allowed_count_prefix),
        Some(&blacklist_mask),
    );

    // Assert
    assert_eq!(sum, 4.0);
    assert_eq!(allowed, 2);
    assert_eq!(blacklisted, 2);
}

#[test]
fn coverage_sum_and_counts_scans_the_mask_when_allowed_count_prefix_is_missing() {
    // Arrange
    // Per-base coverage is still [1, 2, 3, 4]
    // Query [1, 4) -> values [2, 3, 4]
    // Only the middle base is allowed, so:
    //   allowed coverage sum = 3
    //   allowed count = 1
    //   blacklisted count = 2
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];
    let psum_allowed = [0.0, 1.0, 1.0, 4.0, 4.0];
    let blacklist_mask = [0_u8, 1, 0, 1];

    // Act
    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
        1,
        4,
        true,
        &psum_all,
        Some(&psum_allowed),
        None,
        Some(&blacklist_mask),
    );

    // Assert
    assert_eq!(sum, 3.0);
    assert_eq!(allowed, 1);
    assert_eq!(blacklisted, 2);
}

#[test]
fn coverage_sum_and_counts_falls_back_to_full_span_when_mask_support_is_missing() {
    // Arrange
    // Without unmasked prefixes or a blacklist mask, masked mode cannot subtract anything.
    // The helper therefore falls back to treating the whole span as allowed.
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];

    // Act
    let (sum, allowed, blacklisted) =
        coverage_sum_and_counts(1, 3, true, &psum_all, None, None, None);

    // Assert
    assert_eq!(sum, 5.0);
    assert_eq!(allowed, 2);
    assert_eq!(blacklisted, 0);
}

#[test]
fn coverage_sum_and_counts_uses_full_sum_and_span_when_unmasked() {
    // Arrange
    // Per-base coverage is [1, 2, 3, 4]
    // Query [1, 4) -> values [2, 3, 4]
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];
    let psum_allowed = [0.0, 1.0, 1.0, 4.0, 4.0];
    let allowed_count_prefix = [0_u32, 1, 1, 2, 2];
    let blacklist_mask = [0_u8, 1, 0, 1];

    // Act
    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
        1,
        4,
        false,
        &psum_all,
        Some(&psum_allowed),
        Some(&allowed_count_prefix),
        Some(&blacklist_mask),
    );

    // Assert
    assert_eq!(sum, 9.0);
    assert_eq!(allowed, 3);
    assert_eq!(blacklisted, 0);
}

#[test]
fn finalize_value_returns_zero_for_masked_average_with_no_allowed_positions() {
    // Arrange / Act / Assert
    let value = finalize_value(7.5, 0, 100, true, &CoverageWindowAction::Average);
    assert_eq!(value, 0.0);
}

#[test]
fn finalize_value_returns_zero_for_unmasked_average_with_zero_span() {
    // Arrange / Act / Assert
    let value = finalize_value(7.5, 5, 0, false, &CoverageWindowAction::Average);
    assert_eq!(value, 0.0);
}

#[test]
fn finalize_value_returns_sum_for_total_modes_even_when_denominators_are_zero() {
    // Arrange / Act
    let total = finalize_value(7.5, 0, 0, false, &CoverageWindowAction::Total);
    let grouped_total =
        finalize_value(7.5, 0, 0, true, &CoverageWindowAction::TotalOnUniqueBases);

    // Assert
    assert_eq!(total, 7.5);
    assert_eq!(grouped_total, 7.5);
}

#[test]
fn merge_scaled_positional_tiles_orders_tiles_and_scales_values() -> Result<()> {
    // Arrange
    let temp_dir = TempDir::new()?;
    let out_dir = TempDir::new()?;
    // Intentionally write tiles out of index order and across two chromosomes.
    write_zstd_lines(
        &temp_dir,
        "cov.pos.chr2.1.bedgraph.zst",
        &["chr2\t10\t20\t1.25"],
    )?;
    write_zstd_lines(
        &temp_dir,
        "cov.pos.chr1.1.bedgraph.zst",
        &["chr1\t5\t10\t0.5"],
    )?;
    write_zstd_lines(
        &temp_dir,
        "cov.pos.chr1.0.bedgraph.zst",
        &["chr1\t0\t5\t0.3333", "", "chr1\t10\t15\t2"],
    )?;
    write_zstd_lines(
        &temp_dir,
        "cov.pos.chr2.0.bedgraph.zst",
        &["chr2\t0\t10\t0.25"],
    )?;

    // Act
    let merged_path = merge_scaled_positional_tiles(
        temp_dir.path(),
        out_dir.path(),
        &["chr1".to_string(), "chr2".to_string()],
        "cov.pos",
        "merged.bedgraph.zst",
        3.0,
        false,
        2,
        1,
    )?;
    let text = read_text(&merged_path)?;

    // Assert
    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chr1\t0\t5\t1",
            "chr1\t10\t15\t6",
            "chr1\t5\t10\t1.5",
            "chr2\t0\t10\t0.75",
            "chr2\t10\t20\t3.75",
        ]
    );
    Ok(())
}

#[test]
fn merge_scaled_positional_tiles_keeps_indexed_window_column() -> Result<()> {
    // Arrange
    let temp_dir = TempDir::new()?;
    let out_dir = TempDir::new()?;
    write_zstd_lines(
        &temp_dir,
        "cov.pos.chr1.1.tsv.zst",
        &["chr1\t10\t20\t0.5\t7"],
    )?;
    write_zstd_lines(
        &temp_dir,
        "cov.pos.chr1.0.tsv.zst",
        &["chr1\t0\t10\t1.25\t3"],
    )?;

    // Act
    let merged_path = merge_scaled_positional_tiles(
        temp_dir.path(),
        out_dir.path(),
        &["chr1".to_string()],
        "cov.pos",
        "merged.tsv.zst",
        4.0,
        true,
        3,
        1,
    )?;
    let text = read_text(&merged_path)?;

    // Assert
    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec!["chr1\t0\t10\t5\t3", "chr1\t10\t20\t2\t7"]
    );
    Ok(())
}

#[test]
fn merge_scaled_positional_tiles_is_basewise_invariant_to_tile_segmentation() -> Result<()> {
    // Arrange
    let first_temp_dir = TempDir::new()?;
    let second_temp_dir = TempDir::new()?;
    let first_out_dir = TempDir::new()?;
    let second_out_dir = TempDir::new()?;

    // Same underlying pre-scaled signal on chr1:
    // [0,5) -> 0.5, [5,15) -> 1.0, [15,20) -> 0.25
    write_zstd_lines(
        &first_temp_dir,
        "cov.pos.chr1.0.bedgraph.zst",
        &["chr1\t0\t5\t0.5", "chr1\t5\t10\t1"],
    )?;
    write_zstd_lines(
        &first_temp_dir,
        "cov.pos.chr1.1.bedgraph.zst",
        &["chr1\t10\t15\t1", "chr1\t15\t20\t0.25"],
    )?;

    write_zstd_lines(
        &second_temp_dir,
        "cov.pos.chr1.0.bedgraph.zst",
        &["chr1\t0\t3\t0.5", "chr1\t3\t5\t0.5"],
    )?;
    write_zstd_lines(
        &second_temp_dir,
        "cov.pos.chr1.1.bedgraph.zst",
        &["chr1\t5\t8\t1", "chr1\t8\t15\t1"],
    )?;
    write_zstd_lines(
        &second_temp_dir,
        "cov.pos.chr1.2.bedgraph.zst",
        &["chr1\t15\t18\t0.25", "chr1\t18\t20\t0.25"],
    )?;

    // Act
    let first_merged = merge_scaled_positional_tiles(
        first_temp_dir.path(),
        first_out_dir.path(),
        &["chr1".to_string()],
        "cov.pos",
        "merged_a.bedgraph.zst",
        2.0,
        false,
        2,
        1,
    )?;
    let second_merged = merge_scaled_positional_tiles(
        second_temp_dir.path(),
        second_out_dir.path(),
        &["chr1".to_string()],
        "cov.pos",
        "merged_b.bedgraph.zst",
        2.0,
        false,
        2,
        1,
    )?;
    let first_text = read_text(&first_merged)?;
    let second_text = read_text(&second_merged)?;

    // Assert
    let expected = {
        let mut coverage = vec![0.0; 20];
        for position in 0..5 {
            coverage[position] = 1.0;
        }
        for position in 5..15 {
            coverage[position] = 2.0;
        }
        for position in 15..20 {
            coverage[position] = 0.5;
        }
        coverage
    };
    assert_eq!(dense_bedgraph(&first_text, "chr1", 20), expected);
    assert_eq!(dense_bedgraph(&second_text, "chr1", 20), expected);
    Ok(())
}
