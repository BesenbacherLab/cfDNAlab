use super::{TileAggregateTempFiles, reduce_bed_rows, reduce_size_rows};
use crate::shared::{
    bed::Windows,
    interval::{IndexedInterval, Interval},
};
use anyhow::Result;
use fxhash::FxHashMap;
use std::path::Path;
use tempfile::TempDir;

fn write_text(path: &Path, text: &str) -> Result<()> {
    std::fs::write(path, text)?;
    Ok(())
}

#[test]
fn bed_window_reducer_uses_explicit_tile_outputs_and_ignores_matching_decoys() -> Result<()> {
    // Arrange
    // This is a BED-window aggregate reducer test. The temporary rows are not BED files; they are
    // tile-local aggregate contribution rows keyed by the original BED window index.
    //
    // Intended input contract:
    // - the reducer receives the exact partial-row path and optional cross-index path for each tile
    // - unrelated files in the same temporary directory are ignored, even when their names match
    //   the old discovery convention
    //
    // Window orig_idx 7 spans [10,15). It crosses two tile cores, so both tiles list `7` in their
    // cross-index file. The exact additive reduction is:
    //   coverage_sum = 2 + 3 = 5
    //   eligible_positions = 2 + 3 = 5
    //   blacklisted_positions = 0 + 0 = 0
    let temp_dir = TempDir::new()?;
    let tile0_partials = temp_dir.path().join("tile0.window_rows");
    let tile0_cross = temp_dir.path().join("tile0.cross_keys");
    let tile1_partials = temp_dir.path().join("tile1.window_rows");
    let tile1_cross = temp_dir.path().join("tile1.cross_keys");
    write_text(&tile0_partials, "7\t2\t2\t0\n")?;
    write_text(&tile0_cross, "7\n")?;
    write_text(&tile1_partials, "7\t3\t3\t0\n")?;
    write_text(&tile1_cross, "7\n")?;

    // This decoy is intentionally named like an old discoverable partial file. If the reducer
    // scans the temp directory instead of using the typed paths above, the expected row changes.
    let discoverable_decoy = temp_dir.path().join("run.part.chrom-000000.9.tsv");
    write_text(&discoverable_decoy, "7\t999\t999\t0\n")?;

    let tile_outputs = vec![
        TileAggregateTempFiles {
            tile_index: 0,
            partials_path: tile0_partials,
            cross_index_path: Some(tile0_cross),
        },
        TileAggregateTempFiles {
            tile_index: 1,
            partials_path: tile1_partials,
            cross_index_path: Some(tile1_cross),
        },
    ];
    let windows = vec![IndexedInterval::new(10_u64, 15_u64, 7_u64)?];
    let mut reduced_rows = Vec::new();

    // Act
    reduce_bed_rows("chr1", &tile_outputs, &windows, false, |row| {
        reduced_rows.push(row);
        Ok(())
    })?;

    // Assert
    assert_eq!(reduced_rows.len(), 1);
    let row = reduced_rows[0];
    assert_eq!(row.idx, 7);
    assert_eq!(row.interval, Interval::new(10_u64, 15_u64)?);
    assert_eq!(row.coverage_sum, 5.0);
    assert_eq!(row.eligible_positions, 5);
    assert_eq!(row.blacklisted_positions, 0);
    assert_eq!(row.nonzero_positions, 0);
    assert_eq!(row.coverage_sum_of_squares, 0.0);

    Ok(())
}

#[test]
fn grouped_segment_reducer_uses_explicit_tile_outputs_before_group_folding() -> Result<()> {
    // Arrange
    // Grouped aggregate output first reduces internal segment rows, then folds those reduced
    // segment rows into group rows. This test pins the first stage only: the reducer must consume
    // the explicit tile output paths for grouped segment partial rows, not discover files by name.
    //
    // Segment idx 3 belongs to a grouped layout elsewhere. Its span is [20,30), and it crosses two
    // tile cores. The exact additive reduction is:
    //   coverage_sum = 4 + 6 = 10
    //   eligible_positions = 4 + 6 = 10
    //   blacklisted_positions = 0
    let temp_dir = TempDir::new()?;
    let tile0_partials = temp_dir.path().join("tile0.group_segment_rows");
    let tile0_cross = temp_dir.path().join("tile0.group_segment_cross_keys");
    let tile1_partials = temp_dir.path().join("tile1.group_segment_rows");
    let tile1_cross = temp_dir.path().join("tile1.group_segment_cross_keys");
    write_text(&tile0_partials, "3\t4\t4\t0\n")?;
    write_text(&tile0_cross, "3\n")?;
    write_text(&tile1_partials, "3\t6\t6\t0\n")?;
    write_text(&tile1_cross, "3\n")?;

    let discoverable_decoy = temp_dir.path().join("run.part.chrom-000000.9.tsv");
    write_text(&discoverable_decoy, "3\t999\t999\t0\n")?;

    let tile_outputs = vec![
        TileAggregateTempFiles {
            tile_index: 0,
            partials_path: tile0_partials,
            cross_index_path: Some(tile0_cross),
        },
        TileAggregateTempFiles {
            tile_index: 1,
            partials_path: tile1_partials,
            cross_index_path: Some(tile1_cross),
        },
    ];
    let grouped_segments = vec![IndexedInterval::new(20_u64, 30_u64, 3_u64)?];
    let mut reduced_rows = Vec::new();

    // Act
    reduce_bed_rows("chr1", &tile_outputs, &grouped_segments, false, |row| {
        reduced_rows.push(row);
        Ok(())
    })?;

    // Assert
    assert_eq!(reduced_rows.len(), 1);
    let row = reduced_rows[0];
    assert_eq!(row.idx, 3);
    assert_eq!(row.interval, Interval::new(20_u64, 30_u64)?);
    assert_eq!(row.coverage_sum, 10.0);
    assert_eq!(row.eligible_positions, 10);
    assert_eq!(row.blacklisted_positions, 0);
    assert_eq!(row.nonzero_positions, 0);
    assert_eq!(row.coverage_sum_of_squares, 0.0);

    Ok(())
}

#[test]
fn bed_summary_reducer_uses_explicit_tile_outputs_for_raw_moment_rows() -> Result<()> {
    // Arrange
    // Summary-stat aggregate contribution rows carry the same BED-window identity as basic rows
    // plus raw-moment columns. The explicit-path contract must apply to this schema too.
    //
    // Window orig_idx 11 spans [100,110) and receives two tile contributions:
    //   coverage_sum = 1.5 + 2.5 = 4.0
    //   eligible_positions = 4 + 6 = 10
    //   nonzero_positions = 2 + 3 = 5
    //   coverage_sum_of_squares = 2.25 + 6.25 = 8.5
    let temp_dir = TempDir::new()?;
    let tile0_partials = temp_dir.path().join("tile0.summary_window_rows");
    let tile0_cross = temp_dir.path().join("tile0.summary_cross_keys");
    let tile1_partials = temp_dir.path().join("tile1.summary_window_rows");
    let tile1_cross = temp_dir.path().join("tile1.summary_cross_keys");
    write_text(&tile0_partials, "11\t1.5\t4\t0\t2\t2.25\n")?;
    write_text(&tile0_cross, "11\n")?;
    write_text(&tile1_partials, "11\t2.5\t6\t0\t3\t6.25\n")?;
    write_text(&tile1_cross, "11\n")?;
    write_text(
        &temp_dir.path().join("run.part.chrom-000000.9.tsv"),
        "11\t999\t999\t0\t999\t999\n",
    )?;

    let tile_outputs = vec![
        TileAggregateTempFiles {
            tile_index: 0,
            partials_path: tile0_partials,
            cross_index_path: Some(tile0_cross),
        },
        TileAggregateTempFiles {
            tile_index: 1,
            partials_path: tile1_partials,
            cross_index_path: Some(tile1_cross),
        },
    ];
    let windows = vec![IndexedInterval::new(100_u64, 110_u64, 11_u64)?];
    let mut reduced_rows = Vec::new();

    // Act
    reduce_bed_rows("chr1", &tile_outputs, &windows, true, |row| {
        reduced_rows.push(row);
        Ok(())
    })?;

    // Assert
    assert_eq!(reduced_rows.len(), 1);
    let row = reduced_rows[0];
    assert_eq!(row.idx, 11);
    assert_eq!(row.interval, Interval::new(100_u64, 110_u64)?);
    assert_eq!(row.coverage_sum, 4.0);
    assert_eq!(row.eligible_positions, 10);
    assert_eq!(row.blacklisted_positions, 0);
    assert_eq!(row.nonzero_positions, 5);
    assert_eq!(row.coverage_sum_of_squares, 8.5);

    Ok(())
}

#[test]
fn size_reducer_uses_explicit_tile_outputs_and_logical_bin_start() -> Result<()> {
    // Arrange
    // Logical fixed-size bin [40,80) crosses a tile boundary and chromosome end is 75, so the
    // final interval must be clipped to [40,75). The reducer identity is the full logical
    // bin start, not a clipped tile-local span.
    //
    // Hand-derived additive totals:
    //   coverage_sum = 8 + 32 = 40
    //   eligible_positions = 8 + 32 = 40
    //   blacklisted_positions = 0
    let temp_dir = TempDir::new()?;
    let tile0_partials = temp_dir.path().join("tile0.size_rows");
    let tile0_cross = temp_dir.path().join("tile0.size_cross_keys");
    let tile1_partials = temp_dir.path().join("tile1.size_rows");
    let tile1_cross = temp_dir.path().join("tile1.size_cross_keys");
    write_text(&tile0_partials, "40\t80\t8\t8\t0\n")?;
    write_text(&tile0_cross, "40\n")?;
    write_text(&tile1_partials, "40\t80\t32\t32\t0\n")?;
    write_text(&tile1_cross, "40\n")?;

    // Decoy matching the old filename-based discovery convention. Explicit reducer inputs should
    // make this file invisible to the reduction.
    let discoverable_decoy = temp_dir.path().join("run.part.chrom-000000.9.tsv");
    write_text(&discoverable_decoy, "40\t80\t999\t999\t0\n")?;

    let tile_outputs = vec![
        TileAggregateTempFiles {
            tile_index: 0,
            partials_path: tile0_partials,
            cross_index_path: Some(tile0_cross),
        },
        TileAggregateTempFiles {
            tile_index: 1,
            partials_path: tile1_partials,
            cross_index_path: Some(tile1_cross),
        },
    ];
    let mut reduced_rows = Vec::new();

    // Act
    reduce_size_rows(
        "chr1",
        &tile_outputs,
        75,
        false,
        |row| {
            reduced_rows.push(row);
            Ok(())
        },
    )?;

    // Assert
    assert_eq!(reduced_rows.len(), 1);
    let row = reduced_rows[0];
    assert_eq!(row.idx, 40);
    assert_eq!(row.interval, Interval::new(40_u64, 75_u64)?);
    assert_eq!(row.coverage_sum, 40.0);
    assert_eq!(row.eligible_positions, 40);
    assert_eq!(row.blacklisted_positions, 0);
    assert_eq!(row.nonzero_positions, 0);
    assert_eq!(row.coverage_sum_of_squares, 0.0);

    Ok(())
}

#[test]
fn size_summary_reducer_uses_explicit_tile_outputs_for_raw_moment_rows() -> Result<()> {
    // Arrange
    // Size summary-stat rows use logical bin start/end plus raw moments. The reducer must use the
    // explicit paths and continue grouping by logical bin start.
    //
    // Bin [80,120) crosses two tiles and chromosome end is 110, so final interval is [80,110).
    // Exact additive totals:
    //   coverage_sum = 3 + 7 = 10
    //   eligible_positions = 10 + 20 = 30
    //   nonzero_positions = 4 + 6 = 10
    //   coverage_sum_of_squares = 5 + 13 = 18
    let temp_dir = TempDir::new()?;
    let tile0_partials = temp_dir.path().join("tile0.summary_size_rows");
    let tile0_cross = temp_dir.path().join("tile0.summary_size_cross_keys");
    let tile1_partials = temp_dir.path().join("tile1.summary_size_rows");
    let tile1_cross = temp_dir.path().join("tile1.summary_size_cross_keys");
    write_text(&tile0_partials, "80\t120\t3\t10\t0\t4\t5\n")?;
    write_text(&tile0_cross, "80\n")?;
    write_text(&tile1_partials, "80\t120\t7\t20\t0\t6\t13\n")?;
    write_text(&tile1_cross, "80\n")?;
    write_text(
        &temp_dir.path().join("run.part.chrom-000000.9.tsv"),
        "80\t120\t999\t999\t0\t999\t999\n",
    )?;

    let tile_outputs = vec![
        TileAggregateTempFiles {
            tile_index: 0,
            partials_path: tile0_partials,
            cross_index_path: Some(tile0_cross),
        },
        TileAggregateTempFiles {
            tile_index: 1,
            partials_path: tile1_partials,
            cross_index_path: Some(tile1_cross),
        },
    ];
    let mut reduced_rows = Vec::new();

    // Act
    reduce_size_rows(
        "chr1",
        &tile_outputs,
        110,
        true,
        |row| {
            reduced_rows.push(row);
            Ok(())
        },
    )?;

    // Assert
    assert_eq!(reduced_rows.len(), 1);
    let row = reduced_rows[0];
    assert_eq!(row.idx, 80);
    assert_eq!(row.interval, Interval::new(80_u64, 110_u64)?);
    assert_eq!(row.coverage_sum, 10.0);
    assert_eq!(row.eligible_positions, 30);
    assert_eq!(row.blacklisted_positions, 0);
    assert_eq!(row.nonzero_positions, 10);
    assert_eq!(row.coverage_sum_of_squares, 18.0);

    Ok(())
}

#[test]
fn size_reducer_accepts_missing_cross_index_for_single_contribution_partials() -> Result<()> {
    // Arrange
    // Aligned restore-mean by-size tiles are intended to write raw partial rows without a
    // cross-index file. Missing cross-index still means one expected contribution for each
    // logical bin.
    let temp_dir = TempDir::new()?;
    let tile0_partials = temp_dir.path().join("tile0.aligned_size_rows");
    write_text(&tile0_partials, "0\t40\t20\t40\t0\n")?;

    let tile_outputs = vec![TileAggregateTempFiles {
        tile_index: 0,
        partials_path: tile0_partials,
        cross_index_path: None,
    }];
    let mut reduced_rows = Vec::new();

    // Act
    reduce_size_rows(
        "chr1",
        &tile_outputs,
        40,
        false,
        |row| {
            reduced_rows.push(row);
            Ok(())
        },
    )?;

    // Assert
    assert_eq!(reduced_rows.len(), 1);
    let row = reduced_rows[0];
    assert_eq!(row.idx, 0);
    assert_eq!(row.interval, Interval::new(0_u64, 40_u64)?);
    assert_eq!(row.coverage_sum, 20.0);
    assert_eq!(row.eligible_positions, 40);
    assert_eq!(row.blacklisted_positions, 0);
    assert_eq!(row.nonzero_positions, 0);
    assert_eq!(row.coverage_sum_of_squares, 0.0);

    Ok(())
}

#[test]
fn reducer_rejects_duplicate_explicit_tile_indices() -> Result<()> {
    // Arrange
    // The explicit tile-output list should be a typed replacement for one tile's outputs, not a
    // loose bag of paths. Duplicate tile indices make diagnostics and contribution counts
    // ambiguous, so they should fail before any reduction is attempted.
    let temp_dir = TempDir::new()?;
    let first_partials = temp_dir.path().join("first.rows");
    let second_partials = temp_dir.path().join("second.rows");
    write_text(&first_partials, "0\t1\t1\t0\n")?;
    write_text(&second_partials, "0\t1\t1\t0\n")?;

    let tile_outputs = vec![
        TileAggregateTempFiles {
            tile_index: 0,
            partials_path: first_partials,
            cross_index_path: None,
        },
        TileAggregateTempFiles {
            tile_index: 0,
            partials_path: second_partials,
            cross_index_path: None,
        },
    ];
    let windows = vec![IndexedInterval::new(0_u64, 1_u64, 0_u64)?];

    // Act
    let err = reduce_bed_rows("chr1", &tile_outputs, &windows, false, |_| Ok(()))
        .expect_err("duplicate tile indices should fail");

    // Assert
    assert!(
        err.to_string().contains("duplicate tile index 0"),
        "unexpected error: {err:#}"
    );

    Ok(())
}

#[test]
fn reducer_reports_missing_explicit_partial_path() -> Result<()> {
    // Arrange
    // A missing returned partial path is a broken tile-output contract. The reducer should report
    // the concrete path it was asked to read instead of silently falling back to discovery.
    let temp_dir = TempDir::new()?;
    let missing_partials = temp_dir.path().join("missing_tile_rows");
    let tile_outputs = vec![TileAggregateTempFiles {
        tile_index: 0,
        partials_path: missing_partials.clone(),
        cross_index_path: None,
    }];
    let windows = vec![IndexedInterval::new(0_u64, 1_u64, 0_u64)?];

    // Act
    let err = reduce_bed_rows("chr1", &tile_outputs, &windows, false, |_| Ok(()))
        .expect_err("missing explicit partial path should fail");

    // Assert
    let message = err.to_string();
    assert!(
        message.contains("missing_tile_rows") || message.contains(&missing_partials.display().to_string()),
        "unexpected error: {err:#}"
    );

    Ok(())
}

#[test]
fn bed_reducer_partitions_explicit_tile_outputs_by_requested_chromosome_order() -> Result<()> {
    // Arrange
    // Multi-chromosome aggregate reduction should use the caller-provided chromosome grouping and
    // requested chromosome order. It must not infer chromosome identity from filenames or temp-dir
    // order.
    let temp_dir = TempDir::new()?;
    let chr1_partials = temp_dir.path().join("first_chromosome_rows");
    let chr2_partials = temp_dir.path().join("second_chromosome_rows");
    write_text(&chr1_partials, "1\t10\t10\t0\n")?;
    write_text(&chr2_partials, "2\t20\t20\t0\n")?;

    let mut tile_outputs_by_chr = FxHashMap::default();
    tile_outputs_by_chr.insert(
        "chr1".to_string(),
        vec![TileAggregateTempFiles {
            tile_index: 0,
            partials_path: chr1_partials,
            cross_index_path: None,
        }],
    );
    tile_outputs_by_chr.insert(
        "chr2".to_string(),
        vec![TileAggregateTempFiles {
            tile_index: 0,
            partials_path: chr2_partials,
            cross_index_path: None,
        }],
    );

    let mut windows_by_chr = FxHashMap::default();
    windows_by_chr.insert("chr1".to_string(), Windows::from_tuples(&[(10, 20, 1)])?);
    windows_by_chr.insert("chr2".to_string(), Windows::from_tuples(&[(30, 50, 2)])?);

    let mut observed = Vec::new();

    // Act
    for chromosome in ["chr2".to_string(), "chr1".to_string()] {
        let Some(windows_for_chr) = windows_by_chr.get(&chromosome) else {
            continue;
        };
        let empty_outputs = Vec::new();
        let tile_outputs = tile_outputs_by_chr
            .get(&chromosome)
            .unwrap_or(&empty_outputs);
        reduce_bed_rows(
            &chromosome,
            tile_outputs,
            windows_for_chr.as_slice(),
            false,
            |row| {
                observed.push((chromosome.clone(), row));
                Ok(())
            },
        )?;
    }

    // Assert
    assert_eq!(observed.len(), 2);
    assert_eq!(observed[0].0, "chr2");
    assert_eq!(observed[0].1.idx, 2);
    assert_eq!(observed[0].1.interval, Interval::new(30_u64, 50_u64)?);
    assert_eq!(observed[0].1.coverage_sum, 20.0);
    assert_eq!(observed[1].0, "chr1");
    assert_eq!(observed[1].1.idx, 1);
    assert_eq!(observed[1].1.interval, Interval::new(10_u64, 20_u64)?);
    assert_eq!(observed[1].1.coverage_sum, 10.0);

    Ok(())
}
