#![cfg(feature = "cmd_lengths")]
//! Public API tests for Rust output loaders for `cfdna lengths`.
//!
//! These tests exercise the public crate boundary: length-count outputs should
//! load without exposing TSV schema or compression details.

use cfdnalab::{
    interval::Interval,
    output_loaders::{LengthOutputMode, LengthRowMetadata, load_lengths_output},
};
use flate2::{Compression, write::GzEncoder};
use std::{fs::File, io::Write, path::Path};
use tempfile::TempDir;

/// Verify that a plain global lengths TSV loads with length-bin metadata.
#[test]
fn load_lengths_output_reads_global_plain_tsv() -> anyhow::Result<()> {
    // Arrange:
    // Global output has no row metadata. The count headers encode one single-bp
    // bin and one wider half-open bin.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(&path, "count_30\tcount_31_40\n12\t3.5\n")?;

    // Act
    let loaded = load_lengths_output(&path)?;

    // Assert
    assert_eq!(loaded.row_mode(), LengthOutputMode::Global);
    assert_eq!(loaded.row_metadata(), &LengthRowMetadata::Global);
    assert_eq!(loaded.row_count(), 1);
    assert_eq!(loaded.length_bin_count(), 2);
    assert_eq!(loaded.counts().shape(), (1, 2));
    assert_eq!(loaded.count(0, 0), Some(12.0));
    assert_eq!(loaded.count(0, 1), Some(3.5));
    assert_eq!(loaded.length_bin_for_length(30), Some(0));
    assert_eq!(loaded.length_bin_for_length(35), Some(1));
    assert_eq!(loaded.length_bin_for_length(40), None);
    assert_eq!(loaded.length_bins()[0].as_tuple(), (30, 31, 0));
    assert_eq!(loaded.length_bins()[1].as_tuple(), (31, 40, 1));
    Ok(())
}

/// Verify half-open overlap selection for fragment length bins.
#[test]
fn lengths_output_selects_length_bins_overlapping_half_open_range() -> anyhow::Result<()> {
    // Arrange:
    // A half-open query range selects whole length bins that overlap it,
    // including partial overlap at either end but excluding bins that only touch
    // the query boundary.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(&path, "count_30_60\tcount_60_90\tcount_90_120\n1\t2\t3\n")?;
    let loaded = load_lengths_output(&path)?;

    // Act
    let broad_selection = loaded.length_bins_overlapping_range(Interval::new(55, 95)?)?;
    let exact_second_bin = loaded.length_bins_overlapping_range(Interval::new(60, 90)?)?;
    let left_edge_selection = loaded.length_bins_overlapping_range(Interval::new(59, 60)?)?;
    let right_edge_selection = loaded.length_bins_overlapping_range(Interval::new(120, 130)?);

    // Assert
    assert_eq!(
        broad_selection
            .iter()
            .map(|bin| bin.as_tuple())
            .collect::<Vec<_>>(),
        vec![(30, 60, 0), (60, 90, 1), (90, 120, 2)]
    );
    assert_eq!(
        exact_second_bin
            .iter()
            .map(|bin| bin.as_tuple())
            .collect::<Vec<_>>(),
        vec![(60, 90, 1)]
    );
    assert_eq!(
        left_edge_selection
            .iter()
            .map(|bin| bin.as_tuple())
            .collect::<Vec<_>>(),
        vec![(30, 60, 0)]
    );
    assert!(
        right_edge_selection
            .expect_err("touching range should not overlap")
            .to_string()
            .contains("does not overlap any length bins")
    );
    Ok(())
}

/// Verify that windowed lengths rows keep blacklist metadata beside counts.
#[test]
fn load_lengths_output_reads_window_rows_with_blacklist_fraction() -> anyhow::Result<()> {
    // Arrange:
    // Windowed output is keyed by `chrom/start/end`. The optional blacklist
    // fraction should stay with the row metadata, not the numeric count matrix.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(
        &path,
        concat!(
            "chrom\tstart\tend\tblacklisted_fraction\tcount_30\tcount_31_40\n",
            "chr1\t0\t100\t0.25\t12\t3.5\n",
            "chr2\t50\t75\t0\t0.25\t7\n",
        ),
    )?;

    // Act
    let loaded = load_lengths_output(&path)?;

    // Assert
    assert_eq!(loaded.row_mode(), LengthOutputMode::Windows);
    assert_eq!(loaded.counts().shape(), (2, 2));
    assert_eq!(loaded.counts_row_major(), &[12.0, 3.5, 0.25, 7.0]);
    assert_eq!(loaded.count(1, 1), Some(7.0));
    let windows = loaded.window_metadata()?;
    assert_eq!(windows[0].index, 0);
    assert_eq!(windows[0].chrom, "chr1");
    assert_eq!(windows[0].interval.as_tuple(), (0, 100));
    assert_eq!(windows[0].blacklisted_fraction, Some(0.25));
    assert_eq!(windows[1].index, 1);
    assert_eq!(windows[1].chrom, "chr2");
    assert_eq!(windows[1].interval.as_tuple(), (50, 75));
    assert_eq!(windows[1].blacklisted_fraction, Some(0.0));
    let second_window = loaded.window(1)?.expect("second window should exist");
    assert_eq!(second_window.chrom, "chr2");
    Ok(())
}

/// Verify that grouped lengths rows load from zstd-compressed TSV output.
#[test]
fn load_lengths_output_reads_group_rows_from_zstd_tsv() -> anyhow::Result<()> {
    // Arrange:
    // Grouped output is keyed by group name and eligible window count. This also
    // checks that zstd decompression happens inside the loader.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv.zst");
    write_zstd_text(
        &path,
        concat!(
            "group_name\teligible_windows\tblacklisted_fraction\tcount_30\tcount_31_40\n",
            "groupA\t2\t0.125\t12\t3.5\n",
            "groupWithoutWindows\t0\t0\t0\t0\n",
        ),
    )?;

    // Act
    let loaded = load_lengths_output(&path)?;

    // Assert
    assert_eq!(loaded.row_mode(), LengthOutputMode::Groups);
    assert_eq!(loaded.counts().shape(), (2, 2));
    assert_eq!(loaded.counts_row_major(), &[12.0, 3.5, 0.0, 0.0]);
    let groups = loaded.group_metadata()?;
    assert_eq!(groups[0].index, 0);
    assert_eq!(groups[0].name, "groupA");
    assert_eq!(groups[0].eligible_windows, 2);
    assert_eq!(groups[0].blacklisted_fraction, Some(0.125));
    assert_eq!(groups[1].index, 1);
    assert_eq!(groups[1].name, "groupWithoutWindows");
    assert_eq!(groups[1].eligible_windows, 0);
    assert_eq!(groups[1].blacklisted_fraction, Some(0.0));
    let first_group = loaded.group(0)?.expect("first group should exist");
    assert_eq!(first_group.name, "groupA");
    assert_eq!(loaded.group_index("groupWithoutWindows")?, 1);
    assert!(loaded.has_group("groupA"));
    assert!(!loaded.has_group("missing"));
    Ok(())
}

/// Verify that gzip windowed lengths output can omit blacklist fractions.
#[test]
fn load_lengths_output_reads_gzip_window_rows_without_blacklist_fraction() -> anyhow::Result<()> {
    // Arrange:
    // The blacklist column is optional for windowed length outputs, and gzip
    // input should be handled by the loader rather than by downstream code.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv.gz");
    write_gzip_text(
        &path,
        concat!(
            "chrom\tstart\tend\tcount_30\tcount_31_40\n",
            "chr1\t0\t10\t1\t2\n",
        ),
    )?;

    // Act
    let loaded = load_lengths_output(&path)?;

    // Assert
    assert_eq!(loaded.row_mode(), LengthOutputMode::Windows);
    assert_eq!(loaded.counts_row_major(), &[1.0, 2.0]);
    let window = loaded.window(0)?.expect("first window should exist");
    assert_eq!(window.chrom, "chr1");
    assert_eq!(window.interval.as_tuple(), (0, 10));
    assert_eq!(window.blacklisted_fraction, None);
    Ok(())
}

/// Verify row and length-bin selection for windowed lengths output.
#[test]
fn lengths_output_selects_multiple_window_rows_and_length_bins() -> anyhow::Result<()> {
    // Arrange:
    // The selection API should return an owned dense matrix in the requested
    // row and length-bin order, without requiring manual TSV column lookup.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(
        &path,
        concat!(
            "chrom\tstart\tend\tcount_30\tcount_31_40\tcount_40_50\n",
            "chr1\t0\t100\t1\t2\t3\n",
            "chr1\t100\t200\t4\t5\t6\n",
            "chr2\t0\t100\t7\t8\t9\n",
        ),
    )?;
    let loaded = load_lengths_output(&path)?;

    // Act
    let selected = loaded
        .select()
        .windows(&[2, 0])
        .length_bins(&[1, 2])
        .read()?;

    // Assert
    assert_eq!(selected.row_indices(), &[2, 0]);
    assert_eq!(
        selected
            .window_metadata()?
            .iter()
            .map(|window| (window.chrom.as_str(), window.interval.as_tuple()))
            .collect::<Vec<_>>(),
        vec![("chr2", (0, 100)), ("chr1", (0, 100))]
    );
    assert_eq!(
        selected
            .length_bins()
            .iter()
            .map(|bin| bin.as_tuple())
            .collect::<Vec<_>>(),
        vec![(31, 40, 1), (40, 50, 2)]
    );
    assert_eq!(selected.shape(), (2, 2));
    assert_eq!(selected.row_count(), 2);
    assert_eq!(selected.length_bin_count(), 2);
    assert_eq!(selected.counts_row_major(), &[8.0, 9.0, 2.0, 3.0]);
    assert_eq!(selected.count(1, 0), Some(2.0));
    Ok(())
}

/// Verify grouped lengths selection by public group labels.
#[test]
fn lengths_output_selects_group_counts_by_name() -> anyhow::Result<()> {
    // Arrange:
    // Group-name selection resolves public group labels to count-matrix rows
    // and keeps the requested label order.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(
        &path,
        concat!(
            "group_name\teligible_windows\tcount_30\tcount_31_40\tcount_40_50\n",
            "alpha\t2\t1\t2\t3\n",
            "beta\t3\t4\t5\t6\n",
            "gamma\t1\t7\t8\t9\n",
        ),
    )?;
    let loaded = load_lengths_output(&path)?;

    // Act
    let selected = loaded
        .select()
        .groups_by_name(&["gamma", "alpha"])
        .length_bins(&[2, 0])
        .read()?;

    // Assert
    assert_eq!(selected.row_indices(), &[2, 0]);
    assert_eq!(
        selected
            .group_metadata()?
            .iter()
            .map(|group| group.name.as_str())
            .collect::<Vec<_>>(),
        vec!["gamma", "alpha"]
    );
    assert_eq!(
        selected
            .length_bins()
            .iter()
            .map(|bin| bin.as_tuple())
            .collect::<Vec<_>>(),
        vec![(40, 50, 2), (30, 31, 0)]
    );
    assert_eq!(selected.counts_row_major(), &[9.0, 7.0, 3.0, 1.0]);
    assert_eq!(
        selected
            .counts()
            .rows()
            .map(|selected_length_counts| selected_length_counts.iter().copied().sum::<f64>())
            .collect::<Vec<_>>(),
        vec![16.0, 4.0]
    );
    Ok(())
}

/// Verify default selectors and range selectors for grouped lengths output.
#[test]
fn lengths_output_selects_all_bins_or_length_ranges() -> anyhow::Result<()> {
    // Arrange:
    // Optional selectors should support the common cases directly: selected
    // groups across all length bins, all groups across a length range, and all
    // rows/all bins without spelling out either axis.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(
        &path,
        concat!(
            "group_name\teligible_windows\tcount_30_60\tcount_60_90\tcount_90_120\n",
            "alpha\t2\t1\t2\t3\n",
            "beta\t3\t4\t5\t6\n",
            "gamma\t1\t7\t8\t9\n",
        ),
    )?;
    let loaded = load_lengths_output(&path)?;

    // Act
    let selected_named_groups = loaded.select().groups_by_name(&["gamma", "alpha"]).read()?;
    let selected_range = loaded
        .select()
        .groups(&[0, 1, 2])
        .length_range(Interval::new(55, 95)?)
        .read()?;
    let selected_everything = loaded.select().read()?;

    // Assert
    assert_eq!(selected_named_groups.row_indices(), &[2, 0]);
    assert_eq!(
        selected_named_groups
            .group_metadata()?
            .iter()
            .map(|group| group.name.as_str())
            .collect::<Vec<_>>(),
        vec!["gamma", "alpha"]
    );
    assert_eq!(
        selected_named_groups
            .length_bins()
            .iter()
            .map(|bin| bin.as_tuple())
            .collect::<Vec<_>>(),
        vec![(30, 60, 0), (60, 90, 1), (90, 120, 2)]
    );
    assert_eq!(
        selected_named_groups.counts_row_major(),
        &[7.0, 8.0, 9.0, 1.0, 2.0, 3.0]
    );
    assert_eq!(selected_range.row_indices(), &[0, 1, 2]);
    assert_eq!(
        selected_range
            .length_bins()
            .iter()
            .map(|bin| bin.as_tuple())
            .collect::<Vec<_>>(),
        vec![(30, 60, 0), (60, 90, 1), (90, 120, 2)]
    );
    assert_eq!(
        selected_range.counts_row_major(),
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    );
    assert_eq!(selected_everything.row_indices(), &[0, 1, 2]);
    assert_eq!(
        selected_everything.counts_row_major(),
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    );
    Ok(())
}

/// Verify selection errors for row modes, bounds, duplicates, and conflicts.
#[test]
fn lengths_output_selection_reports_wrong_mode_and_bad_indices() -> anyhow::Result<()> {
    // Arrange:
    // Mode-specific selectors should fail clearly instead of returning empty or
    // silently misinterpreted selections.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(&path, "count_30\tcount_31_40\n12\t3.5\n")?;
    let loaded = load_lengths_output(&path)?;

    // Act
    let window_error = loaded
        .select()
        .windows(&[0])
        .length_bins(&[0])
        .read()
        .expect_err("global output is not windowed");
    let row_error = loaded
        .select()
        .rows(&[1])
        .length_bins(&[0])
        .read()
        .expect_err("row index should be validated");
    let length_bin_error = loaded
        .select()
        .rows(&[0])
        .length_bins(&[2])
        .read()
        .expect_err("length bin index should be validated");
    let duplicate_row_error = loaded
        .select()
        .rows(&[0, 0])
        .length_bins(&[0])
        .read()
        .expect_err("duplicate row indices should be rejected");
    let duplicate_bin_error = loaded
        .select()
        .rows(&[0])
        .length_bins(&[1, 1])
        .read()
        .expect_err("duplicate length bin indices should be rejected");
    let conflicting_row_selector_error = loaded
        .select()
        .rows(&[0])
        .windows(&[0])
        .read()
        .expect_err("conflicting row selectors should fail");
    let conflicting_length_selector_error = loaded
        .select()
        .length_bins(&[0])
        .length_range(Interval::new(30, 40)?)
        .read()
        .expect_err("conflicting length selectors should fail");

    // Assert
    assert!(window_error.to_string().contains("not windowed"));
    assert!(row_error.to_string().contains("row index 1 is outside"));
    assert!(
        length_bin_error
            .to_string()
            .contains("length bin index 2 is outside")
    );
    assert!(
        duplicate_row_error
            .to_string()
            .contains("row indices contain duplicate value 0")
    );
    assert!(
        duplicate_bin_error
            .to_string()
            .contains("length bin indices contain duplicate value 1")
    );
    assert!(
        conflicting_row_selector_error
            .to_string()
            .contains("cannot combine rows() and windows() on the row axis")
    );
    assert!(
        conflicting_length_selector_error.to_string().contains(
            "cannot combine length_bins() and length_range() on the fragment length axis"
        )
    );
    Ok(())
}

/// Verify group-name lookup errors for missing and duplicate group labels.
#[test]
fn lengths_output_group_name_lookup_reports_missing_and_duplicate_names() -> anyhow::Result<()> {
    // Arrange:
    // Group labels are user-facing selectors. Missing labels and duplicated
    // labels must be reported rather than resolving to an arbitrary row.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.length_counts.tsv");
    write_text(
        &path,
        concat!(
            "group_name\teligible_windows\tcount_30\n",
            "alpha\t2\t1\n",
            "alpha\t3\t2\n",
            "beta\t1\t3\n",
        ),
    )?;
    let loaded = load_lengths_output(&path)?;

    // Act
    let duplicate_name_error = loaded
        .group_index("alpha")
        .expect_err("duplicate group names should fail");
    let missing_name_error = loaded
        .select()
        .groups_by_name(&["missing"])
        .read()
        .expect_err("missing group name should fail");

    // Assert
    assert!(
        duplicate_name_error
            .to_string()
            .contains("multiple groups named 'alpha'")
    );
    assert!(
        missing_name_error
            .to_string()
            .contains("no group named 'missing'")
    );
    Ok(())
}

/// Verify malformed length-bin headers are rejected.
#[test]
fn load_lengths_output_rejects_malformed_count_headers() -> anyhow::Result<()> {
    // Arrange:
    // Count columns are the length axis. Malformed, gapped, or out-of-bounds
    // headers would make the matrix disagree with the command's length axis.
    let temp = TempDir::new()?;
    let reversed_path = temp.path().join("reversed.length_counts.tsv");
    let gapped_path = temp.path().join("gapped.length_counts.tsv");
    let below_min_path = temp.path().join("below_min.length_counts.tsv");
    let above_max_path = temp.path().join("above_max.length_counts.tsv");
    write_text(&reversed_path, "count_40_31\n1\n")?;
    write_text(&gapped_path, "count_30_40\tcount_50_60\n1\t2\n")?;
    write_text(&below_min_path, "count_9_10\n1\n")?;
    write_text(&above_max_path, "count_50001\n1\n")?;

    // Act
    let reversed_error =
        load_lengths_output(&reversed_path).expect_err("reversed length bin should fail");
    let gapped_error =
        load_lengths_output(&gapped_path).expect_err("gapped length axis should fail");
    let below_min_error =
        load_lengths_output(&below_min_path).expect_err("length bin below minimum should fail");
    let above_max_error =
        load_lengths_output(&above_max_path).expect_err("length bin above maximum should fail");

    // Assert
    assert!(reversed_error.to_string().contains("has end <= start"));
    assert!(
        gapped_error
            .to_string()
            .contains("non-contiguous length count columns")
    );
    assert!(
        below_min_error
            .to_string()
            .contains("starts below minimum supported fragment length")
    );
    assert!(
        above_max_error
            .to_string()
            .contains("ends above maximum supported fragment length edge")
    );
    Ok(())
}

/// Verify length-count values are finite and non-negative.
#[test]
fn load_lengths_output_rejects_invalid_count_values() -> anyhow::Result<()> {
    // Arrange:
    // Length counts come from non-negative fragment weights. The loader should
    // reject corrupt numeric values instead of passing them into downstream
    // sums or plots.
    let temp = TempDir::new()?;
    let negative_count_path = temp.path().join("negative.length_counts.tsv");
    let non_finite_count_path = temp.path().join("non_finite.length_counts.tsv");
    write_text(&negative_count_path, "count_30\n-1\n")?;
    write_text(&non_finite_count_path, "count_30\nNaN\n")?;

    // Act
    let negative_count_error =
        load_lengths_output(&negative_count_path).expect_err("negative count should fail");
    let non_finite_count_error =
        load_lengths_output(&non_finite_count_path).expect_err("non-finite count should fail");

    // Assert
    assert!(
        negative_count_error
            .to_string()
            .contains("outside finite and non-negative range")
    );
    assert!(
        non_finite_count_error
            .to_string()
            .contains("outside finite and non-negative range")
    );
    Ok(())
}

/// Verify missing data rows and malformed length rows are rejected.
#[test]
fn load_lengths_output_rejects_missing_rows_and_malformed_rows() -> anyhow::Result<()> {
    // Arrange:
    // Header-only files, invalid intervals, invalid fractions, and wrong column
    // counts all make the count matrix ambiguous or invalid.
    let temp = TempDir::new()?;
    let header_only_path = temp.path().join("header_only.length_counts.tsv");
    let invalid_interval_path = temp.path().join("invalid_interval.length_counts.tsv");
    let invalid_fraction_path = temp.path().join("invalid_fraction.length_counts.tsv");
    let wrong_column_count_path = temp.path().join("wrong_columns.length_counts.tsv");
    write_text(&header_only_path, "count_30\n")?;
    write_text(
        &invalid_interval_path,
        "chrom\tstart\tend\tcount_30\nchr1\t10\t10\t1\n",
    )?;
    write_text(
        &invalid_fraction_path,
        "group_name\teligible_windows\tblacklisted_fraction\tcount_30\nalpha\t1\t1.5\t1\n",
    )?;
    write_text(
        &wrong_column_count_path,
        "group_name\teligible_windows\tcount_30\nalpha\t1\n",
    )?;

    // Act
    let header_only_error =
        load_lengths_output(&header_only_path).expect_err("header-only file should fail");
    let invalid_interval_error =
        load_lengths_output(&invalid_interval_path).expect_err("invalid interval should fail");
    let invalid_fraction_error =
        load_lengths_output(&invalid_fraction_path).expect_err("invalid fraction should fail");
    let wrong_column_count_error =
        load_lengths_output(&wrong_column_count_path).expect_err("wrong column count should fail");

    // Assert
    assert!(header_only_error.to_string().contains("has no data rows"));
    assert!(
        invalid_interval_error
            .to_string()
            .contains("invalid window interval")
    );
    assert!(
        invalid_fraction_error
            .to_string()
            .contains("outside [0, 1]")
    );
    assert!(
        wrong_column_count_error
            .to_string()
            .contains("has 2 columns, expected 3")
    );
    Ok(())
}

/// Verify bedGraph-like positional input is not accepted as lengths output.
#[test]
fn load_lengths_output_rejects_bedgraph_like_input() -> anyhow::Result<()> {
    // Arrange:
    // A positional track is not a length-count table and should not be guessed
    // into a generic string-column representation.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.per_position.bedgraph");
    write_text(&path, "chr1\t0\t10\t1.5\n")?;

    // Act
    let error = load_lengths_output(&path).expect_err("bedGraph-like input should fail");

    // Assert
    assert!(error.to_string().contains("unsupported header"));
    Ok(())
}

/// Write an uncompressed text fixture.
fn write_text(path: &Path, text: &str) -> anyhow::Result<()> {
    std::fs::write(path, text)?;
    Ok(())
}

/// Write a zstd-compressed text fixture.
fn write_zstd_text(path: &Path, text: &str) -> anyhow::Result<()> {
    let file = File::create(path)?;
    let mut encoder = zstd::Encoder::new(file, 3)?;
    encoder.write_all(text.as_bytes())?;
    encoder.finish()?;
    Ok(())
}

/// Write a gzip-compressed text fixture.
fn write_gzip_text(path: &Path, text: &str) -> anyhow::Result<()> {
    let file = File::create(path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(text.as_bytes())?;
    encoder.finish()?;
    Ok(())
}
