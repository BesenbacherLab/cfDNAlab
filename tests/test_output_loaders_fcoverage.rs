#![cfg(feature = "cmd_fcoverage")]
//! Public API tests for Rust output loaders for non-positional `cfdna fcoverage`.
//!
//! These tests cover aggregate TSV loading only. Positional bedGraph and
//! per-window positional outputs are intentionally rejected by this public API.

use cfdnalab::output_loaders::{
    FCoverageCoefficientOfVariation, FCoverageRowMode, FCoverageValueMode, load_fcoverage_output,
    load_fcoverage_output_with_group_index,
};
use std::{fs::File, io::Write, path::Path};
use tempfile::TempDir;

/// Verify windowed average fcoverage TSV loading and row selection.
#[test]
fn load_fcoverage_output_reads_window_average_tsv() -> anyhow::Result<()> {
    // Arrange:
    // Windowed average output has genomic row identity, one value column, and
    // blacklisted positions. The span is the checked interval length.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.average.tsv");
    write_text(
        &path,
        concat!(
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions\n",
            "chr1\t0\t100\t1.5\t10\n",
            "chr2\t50\t75\t0\t0\n",
        ),
    )?;

    // Act
    let loaded = load_fcoverage_output(&path)?;

    // Assert
    assert_eq!(loaded.row_mode(), FCoverageRowMode::Windows);
    assert_eq!(loaded.value_mode()?, FCoverageValueMode::Average);
    assert_eq!(loaded.signal().label(), "coverage");
    assert_eq!(loaded.row_count(), 2);
    assert_eq!(loaded.values()?, &[1.5, 0.0]);
    assert_eq!(loaded.value(1)?, Some(0.0));

    let windows = loaded.window_metadata()?;
    assert_eq!(windows[0].index, 0);
    assert_eq!(windows[0].chrom, "chr1");
    assert_eq!(windows[0].interval.as_tuple(), (0, 100));
    assert_eq!(windows[0].blacklisted_positions, 10);
    assert_eq!(windows[0].eligible_positions, None);
    assert_eq!(windows[0].blacklisted_fraction(), 0.1);
    let second_window = loaded.window(1)?.expect("second window should exist");
    assert_eq!(second_window.chrom, "chr2");
    assert!(loaded.group_metadata().is_err());

    let selected_rows = loaded.select().windows(&[1, 0]).read()?;
    assert_eq!(selected_rows.row_indices(), &[1, 0]);
    assert_eq!(selected_rows.row_count(), 2);
    assert_eq!(selected_rows.value_mode()?, FCoverageValueMode::Average);
    assert_eq!(selected_rows.signal().label(), "coverage");
    assert_eq!(selected_rows.values()?, &[0.0, 1.5]);
    assert_eq!(selected_rows.value(1)?, Some(1.5));
    assert_eq!(
        selected_rows
            .window_metadata()?
            .iter()
            .map(|window| (window.chrom.as_str(), window.interval.as_tuple()))
            .collect::<Vec<_>>(),
        vec![("chr2", (50, 75)), ("chr1", (0, 100))]
    );
    assert!(selected_rows.group_metadata().is_err());
    Ok(())
}

/// Verify grouped total fcoverage loading from zstd-compressed TSV output.
#[test]
fn load_fcoverage_output_reads_grouped_total_fragment_mass_from_zstd() -> anyhow::Result<()> {
    // Arrange:
    // Length-normalized aggregate outputs use `fragment_mass` in the value
    // header. Grouped rows are keyed by the written group index.
    let temp = TempDir::new()?;
    let path = temp
        .path()
        .join("sample.fcoverage.total_on_unique_bases.tsv.zst");
    let group_index_path = temp.path().join("sample.group_index.tsv");
    write_zstd_text(
        &path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\ttotal_fragment_mass\n",
            "0\t100\t5\t95\t42.5\n",
            "3\t0\t0\t0\t0\n",
        ),
    )?;
    write_text(
        &group_index_path,
        "group_idx\tgroup_name\n0\talpha\n3\tbeta\n",
    )?;

    // Act
    let loaded_without_names = load_fcoverage_output(&path)?;
    let loaded = load_fcoverage_output_with_group_index(&path, &group_index_path)?;

    // Assert
    assert!(loaded_without_names.group_metadata()?[0].name.is_none());
    assert_eq!(loaded.row_mode(), FCoverageRowMode::Groups);
    assert_eq!(loaded.value_mode()?, FCoverageValueMode::Total);
    assert_eq!(loaded.signal().label(), "fragment_mass");
    assert_eq!(loaded.values()?, &[42.5, 0.0]);
    assert!(loaded.window_metadata().is_err());

    let groups = loaded.group_metadata()?;
    assert_eq!(groups[0].index, 0);
    assert_eq!(groups[0].group_idx, 0);
    assert_eq!(groups[0].name.as_deref(), Some("alpha"));
    assert_eq!(groups[0].span_positions, 100);
    assert_eq!(groups[0].blacklisted_positions, 5);
    assert_eq!(groups[0].eligible_positions, 95);
    assert_eq!(groups[0].blacklisted_fraction(), 0.05);
    let first_group = loaded.group(0)?.expect("first group should exist");
    assert_eq!(first_group.group_idx, 0);
    assert_eq!(loaded.group_index("beta")?, 1);
    assert!(loaded.has_group("alpha"));
    assert!(!loaded.has_group("gamma"));
    assert!(groups[1].blacklisted_fraction().is_nan());

    let selected_rows = loaded.select().groups(&[1, 0]).read()?;
    assert_eq!(selected_rows.row_indices(), &[1, 0]);
    assert_eq!(selected_rows.row_count(), 2);
    assert_eq!(selected_rows.value_mode()?, FCoverageValueMode::Total);
    assert_eq!(selected_rows.signal().label(), "fragment_mass");
    assert_eq!(selected_rows.values()?, &[0.0, 42.5]);
    assert_eq!(
        selected_rows
            .group_metadata()?
            .iter()
            .map(|group| group.group_idx)
            .collect::<Vec<_>>(),
        vec![3, 0]
    );
    let selected_rows_by_name = loaded.select().groups_by_name(&["beta", "alpha"]).read()?;
    assert_eq!(selected_rows_by_name.row_indices(), &[1, 0]);
    assert_eq!(selected_rows_by_name.values()?, &[0.0, 42.5]);
    assert_eq!(
        selected_rows_by_name
            .group_metadata()?
            .iter()
            .map(|group| group.name.as_deref())
            .collect::<Vec<_>>(),
        vec![Some("beta"), Some("alpha")]
    );
    Ok(())
}

/// Verify NaN scalar aggregate values are preserved.
#[test]
fn load_fcoverage_output_reads_nan_scalar_aggregate_values() -> anyhow::Result<()> {
    // Arrange:
    // Average fcoverage outputs use NaN when no denominator exists. The loader
    // should preserve that value rather than rejecting it or converting it.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.average.tsv");
    write_text(
        &path,
        concat!(
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions\n",
            "chr1\t0\t40\tNaN\t40\n",
        ),
    )?;

    // Act
    let loaded = load_fcoverage_output(&path)?;
    let selected_rows = loaded.select().windows(&[0]).read()?;

    // Assert
    assert!(loaded.value(0)?.expect("first aggregate value").is_nan());
    assert!(loaded.values()?[0].is_nan());
    assert!(
        selected_rows.values()?[0].is_nan(),
        "selected aggregate NaN should be preserved"
    );
    Ok(())
}

/// Verify window summary stats and thresholded coefficient of variation parsing.
#[test]
fn load_fcoverage_output_reads_window_summary_stats_with_cv_threshold() -> anyhow::Result<()> {
    // Arrange:
    // Summary-stat headers encode the signal label on every metric. The CV
    // column may use the writer's threshold display form.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.summary_stats.tsv");
    write_text(
        &path,
        concat!(
            "chromosome\tstart\tend\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_coverage\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "chr1\t10\t20\t10\t2\t8\t4\t0.5\t12\t20\t1.5\t0.25\t0.5\t>1e6\n",
        ),
    )?;

    // Act
    let loaded = load_fcoverage_output(&path)?;

    // Assert
    assert_eq!(loaded.row_mode(), FCoverageRowMode::Windows);
    assert_eq!(loaded.signal().label(), "coverage");
    assert!(loaded.values().is_err());
    assert!(loaded.value_mode().is_err());

    let windows = loaded.window_metadata()?;
    assert_eq!(windows[0].interval.as_tuple(), (10, 20));
    assert_eq!(windows[0].eligible_positions, Some(8));

    let stats = loaded.summary_stat(0)?.expect("summary row");
    assert_eq!(stats.nonzero_positions, 4);
    assert_eq!(stats.covered_fraction, 0.5);
    assert_eq!(stats.total, 12.0);
    assert_eq!(stats.total_squared, 20.0);
    assert_eq!(stats.average, 1.5);
    assert_eq!(stats.variance, 0.25);
    assert_eq!(stats.sd, 0.5);
    assert_eq!(
        stats.coefficient_of_variation,
        FCoverageCoefficientOfVariation::GreaterThan(1.0e6)
    );

    let selected_rows = loaded.select().windows(&[0]).read()?;
    assert_eq!(selected_rows.row_indices(), &[0]);
    assert_eq!(selected_rows.row_count(), 1);
    assert_eq!(selected_rows.signal().label(), "coverage");
    assert_eq!(selected_rows.summary_stat(0)?, Some(stats));
    assert_eq!(
        selected_rows.window_metadata()?[0].interval.as_tuple(),
        (10, 20)
    );
    Ok(())
}

/// Verify grouped summary-stat selections keep row metadata order.
#[test]
fn fcoverage_grouped_summary_selection_keeps_group_metadata_order() -> anyhow::Result<()> {
    // Arrange:
    // Grouped summary-stat selections should carry selected group metadata in
    // the same order as selected summary rows, so downstream code can zip the
    // two directly.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.summary_stats.tsv");
    write_text(
        &path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_coverage\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "0\t10\t0\t10\t5\t0.5\t15\t25\t1.5\t0.25\t0.5\t0.333\n",
            "3\t20\t2\t18\t9\t0.5\t45\t125\t2.5\t0.75\t0.866\t0.346\n",
        ),
    )?;
    let loaded = load_fcoverage_output(&path)?;

    // Act
    let selected_rows = loaded.select().groups(&[1, 0]).read()?;

    // Assert
    assert_eq!(selected_rows.row_indices(), &[1, 0]);
    assert_eq!(
        selected_rows
            .group_metadata()?
            .iter()
            .map(|group| group.group_idx)
            .collect::<Vec<_>>(),
        vec![3, 0]
    );
    assert_eq!(
        selected_rows
            .summary_stats()?
            .iter()
            .map(|stats| stats.average)
            .collect::<Vec<_>>(),
        vec![2.5, 1.5]
    );
    assert_eq!(
        selected_rows
            .group_metadata()?
            .iter()
            .zip(selected_rows.summary_stats()?)
            .map(|(group, stats)| (group.group_idx, stats.average))
            .collect::<Vec<_>>(),
        vec![(3, 2.5), (0, 1.5)]
    );
    assert!(selected_rows.window_metadata().is_err());
    Ok(())
}

/// Verify named group selection rejects outputs loaded without group names.
#[test]
fn fcoverage_group_name_selection_requires_group_index_file() -> anyhow::Result<()> {
    // Arrange:
    // The aggregate TSV only stores numeric group_idx values. Name-based
    // lookup is only available after loading the matching group-index file.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.total.tsv");
    write_text(
        &path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\ttotal_coverage\n",
            "0\t10\t0\t10\t1\n",
        ),
    )?;
    let loaded = load_fcoverage_output(&path)?;

    // Act
    let error = loaded
        .select()
        .groups_by_name(&["alpha"])
        .read()
        .expect_err("name selection should require a group-index file");
    let empty_name_selection_error = loaded
        .select()
        .groups_by_name::<&str>(&[])
        .read()
        .expect_err("empty name selection should still require a group-index file");

    // Assert
    assert!(error.to_string().contains("no group names loaded"));
    assert!(
        empty_name_selection_error
            .to_string()
            .contains("no group names loaded")
    );
    Ok(())
}

/// Verify NaN values in summary-stat rows are preserved.
#[test]
fn load_fcoverage_output_reads_nan_summary_stats() -> anyhow::Result<()> {
    // Arrange:
    // Summary-stat outputs use NaN for undefined derived statistics, for
    // example when no positions are eligible or the mean is zero.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.summary_stats.tsv");
    write_text(
        &path,
        concat!(
            "chromosome\tstart\tend\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_coverage\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "chr1\t0\t40\t40\t40\t0\t0\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n",
        ),
    )?;

    // Act
    let loaded = load_fcoverage_output(&path)?;
    let selected_rows = loaded.select().windows(&[0]).read()?;

    // Assert
    let stats = loaded.summary_stat(0)?.expect("summary row");
    assert!(stats.covered_fraction.is_nan());
    assert!(stats.total.is_nan());
    assert!(stats.total_squared.is_nan());
    assert!(stats.average.is_nan());
    assert!(stats.variance.is_nan());
    assert!(stats.sd.is_nan());
    let FCoverageCoefficientOfVariation::Value(coefficient_of_variation) =
        stats.coefficient_of_variation
    else {
        panic!("expected ordinary CV value");
    };
    assert!(coefficient_of_variation.is_nan());
    assert!(
        selected_rows
            .summary_stat(0)?
            .expect("selected summary row")
            .average
            .is_nan()
    );
    Ok(())
}

/// Verify duplicate scalar-value row selections are rejected.
#[test]
fn fcoverage_value_selection_rejects_duplicate_rows() -> anyhow::Result<()> {
    // Arrange:
    // Selection rows are explicit output-row indices. Duplicates would make
    // row-aligned downstream data ambiguous, so the loader reports them.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.average.tsv");
    write_text(
        &path,
        concat!(
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions\n",
            "chr1\t0\t100\t1.5\t10\n",
            "chr2\t50\t75\t0\t0\n",
        ),
    )?;
    let loaded = load_fcoverage_output(&path)?;

    // Act
    let error = loaded
        .select()
        .windows(&[0, 0])
        .read()
        .expect_err("duplicate row selection should fail");

    // Assert
    assert!(error.to_string().contains("duplicate value 0"));
    Ok(())
}

/// Verify row-mode, value-mode, and selector-conflict errors.
#[test]
fn fcoverage_selection_reports_wrong_row_or_value_mode() -> anyhow::Result<()> {
    // Arrange:
    // Typed selectors check that the loaded file has the requested row and
    // value mode before returning data.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.summary_stats.tsv");
    write_text(
        &path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_coverage\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "0\t10\t0\t10\t10\t1\t10\t10\t1\t0\t0\tNaN\n",
        ),
    )?;
    let loaded = load_fcoverage_output(&path)?;

    // Act
    let window_error = loaded
        .select()
        .windows(&[0])
        .read()
        .expect_err("grouped output should not provide window rows");
    let selected_rows = loaded.select().groups(&[0]).read()?;
    let value_error = selected_rows
        .values()
        .expect_err("summary-stat output should not provide scalar values");
    let conflicting_row_selector_error = loaded
        .select()
        .rows(&[0])
        .groups(&[0])
        .read()
        .expect_err("conflicting row selectors should fail");

    // Assert
    assert!(window_error.to_string().contains("not windowed"));
    assert!(value_error.to_string().contains("summary stats"));
    assert!(
        conflicting_row_selector_error
            .to_string()
            .contains("cannot combine rows() and groups() on the row axis")
    );
    Ok(())
}

/// Verify summary-stat row selections reject bad row indices.
#[test]
fn fcoverage_selection_reports_bad_summary_row_indices() -> anyhow::Result<()> {
    // Arrange:
    // Summary-stat selections use the same row selector API as scalar values,
    // so out-of-bounds and duplicate rows should fail there too.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.summary_stats.tsv");
    write_text(
        &path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_coverage\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "0\t10\t0\t10\t10\t1\t10\t10\t1\t0\t0\t0.5\n",
        ),
    )?;
    let loaded = load_fcoverage_output(&path)?;

    // Act
    let out_of_bounds_error = loaded
        .select()
        .groups(&[1])
        .read()
        .expect_err("row index should be validated");
    let duplicate_row_error = loaded
        .select()
        .groups(&[0, 0])
        .read()
        .expect_err("duplicate row index should be rejected");

    // Assert
    assert!(
        out_of_bounds_error
            .to_string()
            .contains("group index 1 is outside")
    );
    assert!(
        duplicate_row_error
            .to_string()
            .contains("duplicate value 0")
    );
    Ok(())
}

/// Verify positional fcoverage paths are rejected by the aggregate loader.
#[test]
fn load_fcoverage_output_rejects_positional_bedgraph_paths() -> anyhow::Result<()> {
    // Arrange:
    // Positional files are deliberately outside this loader API.
    let temp = TempDir::new()?;
    let path = temp
        .path()
        .join("sample.fcoverage.per_position.bedgraph.zst");
    write_zstd_text(&path, "chr1\t0\t10\t1\n")?;

    // Act
    let error = load_fcoverage_output(&path).expect_err("positional fcoverage should fail");

    // Assert
    assert!(error.to_string().contains("positional fcoverage outputs"));
    Ok(())
}

/// Verify aggregate files are not rejected just because their prefix has positional words.
#[test]
fn load_fcoverage_output_accepts_aggregate_prefix_with_positional_words() -> anyhow::Result<()> {
    // Arrange:
    // Only the generated positional output suffixes are rejected up front.
    // Other names should be parsed from their aggregate TSV schema.
    let temp = TempDir::new()?;
    let path = temp
        .path()
        .join("sample.per_position.bedgraph.fcoverage.average.tsv");
    write_text(
        &path,
        concat!(
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions\n",
            "chr1\t0\t10\t2\t0\n",
        ),
    )?;

    // Act
    let loaded = load_fcoverage_output(&path)?;

    // Assert
    assert_eq!(loaded.values()?, &[2.0]);
    Ok(())
}

/// Verify empty fcoverage tables and invalid row metadata are rejected.
#[test]
fn load_fcoverage_output_rejects_missing_rows_and_invalid_row_metadata() -> anyhow::Result<()> {
    // Arrange:
    // The loader should reject empty aggregate tables and impossible row
    // metadata before returning a typed output.
    let temp = TempDir::new()?;
    let header_only_path = temp.path().join("header_only.fcoverage.average.tsv");
    let invalid_blacklist_path = temp.path().join("invalid_blacklist.fcoverage.average.tsv");
    let invalid_span_path = temp.path().join("invalid_span.fcoverage.summary_stats.tsv");
    let invalid_group_eligible_path = temp
        .path()
        .join("invalid_group_eligible.fcoverage.total.tsv");
    write_text(
        &header_only_path,
        "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions\n",
    )?;
    write_text(
        &invalid_blacklist_path,
        concat!(
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions\n",
            "chr1\t0\t10\t1\t11\n",
        ),
    )?;
    write_text(
        &invalid_span_path,
        concat!(
            "chromosome\tstart\tend\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_coverage\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "chr1\t0\t10\t11\t0\t10\t10\t1\t10\t10\t1\t0\t0\t0.5\n",
        ),
    )?;
    write_text(
        &invalid_group_eligible_path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\ttotal_coverage\n",
            "0\t10\t0\t11\t1\n",
        ),
    )?;

    // Act
    let header_only_error =
        load_fcoverage_output(&header_only_path).expect_err("header-only file should fail");
    let invalid_blacklist_error =
        load_fcoverage_output(&invalid_blacklist_path).expect_err("bad blacklist should fail");
    let invalid_span_error =
        load_fcoverage_output(&invalid_span_path).expect_err("bad span should fail");
    let invalid_group_eligible_error = load_fcoverage_output(&invalid_group_eligible_path)
        .expect_err("bad eligible positions should fail");

    // Assert
    assert!(header_only_error.to_string().contains("has no data rows"));
    assert!(
        invalid_blacklist_error
            .to_string()
            .contains("blacklisted_positions 11 greater than span 10")
    );
    assert!(
        invalid_span_error
            .to_string()
            .contains("span_positions 11 but interval length 10")
    );
    assert!(
        invalid_group_eligible_error
            .to_string()
            .contains("eligible_positions 11 greater than span_positions 10")
    );
    Ok(())
}

/// Verify invalid value headers and coefficient-of-variation thresholds fail.
#[test]
fn load_fcoverage_output_rejects_bad_value_header_and_cv_threshold() -> anyhow::Result<()> {
    // Arrange:
    // Value headers must identify average or total values, and thresholded CV
    // fields must use a positive finite threshold.
    let temp = TempDir::new()?;
    let unsupported_value_path = temp.path().join("sample.fcoverage.aggregate.tsv");
    let bad_cv_path = temp.path().join("sample.fcoverage.summary_stats.tsv");
    write_text(
        &unsupported_value_path,
        concat!(
            "chromosome\tstart\tend\tmedian_coverage\tblacklisted_positions\n",
            "chr1\t0\t10\t1\t0\n",
        ),
    )?;
    write_text(
        &bad_cv_path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_coverage\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "0\t10\t0\t10\t10\t1\t10\t10\t1\t0\t0\t>0\n",
        ),
    )?;

    // Act
    let unsupported_value_error =
        load_fcoverage_output(&unsupported_value_path).expect_err("value header should fail");
    let bad_cv_error =
        load_fcoverage_output(&bad_cv_path).expect_err("bad CV threshold should fail");

    // Assert
    assert!(
        unsupported_value_error
            .to_string()
            .contains("unsupported aggregate value column 'median_coverage'")
    );
    assert!(
        bad_cv_error
            .to_string()
            .contains("invalid coefficient_of_variation threshold")
    );
    Ok(())
}

/// Verify mixed summary-stat signal labels are rejected.
#[test]
fn load_fcoverage_output_rejects_inconsistent_summary_signal_headers() -> anyhow::Result<()> {
    // Arrange:
    // Mixed signal labels would make the loaded summary stats ambiguous.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.fcoverage.summary_stats.tsv");
    write_text(
        &path,
        concat!(
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\t",
            "nonzero_positions\tcovered_fraction\ttotal_coverage\ttotal_squared_fragment_mass\t",
            "average_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\n",
            "0\t10\t0\t10\t10\t1\t10\t10\t1\t0\t0\tNaN\n",
        ),
    )?;

    // Act
    let error = load_fcoverage_output(&path).expect_err("mixed signal labels should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("inconsistent summary-stat signal labels")
    );
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
