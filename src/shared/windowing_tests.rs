use super::*;
use crate::commands::cli_common::WindowSpec;
use crate::shared::{
    bam::Contigs,
    bed::{GroupedWindows, Windows},
    interval::IndexedInterval,
};
use std::path::PathBuf;

fn contigs_with_lengths(entries: &[(&str, u32)]) -> Contigs {
    let mut contigs = Contigs {
        contigs: FxHashMap::default(),
    };
    for (name, len) in entries {
        contigs.contigs.insert((*name).to_string(), (0_i32, *len));
    }
    contigs
}

#[test]
fn window_context_original_idx_uses_offsets_for_size_windows() {
    // Arrange
    let spec = WindowSpec::Size(100);
    let ctx = WindowContext {
        spec: &spec,
        windows: None,
        chr_idx_offset: 12,
    };

    // Act / Assert
    assert_eq!(ctx.original_idx(3), 15);
}

#[test]
fn window_context_original_idx_uses_embedded_bed_indices() -> Result<()> {
    // Arrange
    let spec = WindowSpec::Bed(PathBuf::from("windows.bed"));
    let windows = vec![IndexedInterval::new(10, 20, 7_u64)?];
    let ctx = WindowContext {
        spec: &spec,
        windows: Some(&windows),
        chr_idx_offset: 0,
    };

    // Act / Assert
    assert_eq!(ctx.original_idx(0), 7);
    Ok(())
}

#[test]
fn compute_window_offsets_counts_size_windows_across_chromosomes() -> Result<()> {
    // Arrange
    let chromosomes = vec!["chr1".to_string(), "chr2".to_string()];
    let contigs = contigs_with_lengths(&[("chr1", 250), ("chr2", 90)]);

    // Act
    let (total, offsets) =
        compute_window_offsets(&WindowSpec::Size(100), &chromosomes, &contigs, None)?;

    // Assert
    assert_eq!(total, 4);
    assert_eq!(offsets.get("chr1"), Some(&0));
    assert_eq!(offsets.get("chr2"), Some(&3));
    Ok(())
}

#[test]
fn compute_window_offsets_preserves_bed_original_indices() -> Result<()> {
    // Arrange
    let chromosomes = vec!["chr1".to_string(), "chr2".to_string()];
    let contigs = contigs_with_lengths(&[("chr1", 100), ("chr2", 100)]);
    let mut windows_map = FxHashMap::default();
    windows_map.insert("chr1".to_string(), Windows::from_tuples(&[(0, 10, 5_u64)])?);
    windows_map.insert(
        "chr2".to_string(),
        Windows::from_tuples(&[(10, 20, 9_u64)])?,
    );

    // Act
    let (total, offsets) = compute_window_offsets(
        &WindowSpec::Bed(PathBuf::from("windows.bed")),
        &chromosomes,
        &contigs,
        Some(&windows_map),
    )?;

    // Assert
    assert_eq!(total, 2);
    assert_eq!(offsets.get("chr1"), Some(&0));
    assert_eq!(offsets.get("chr2"), Some(&0));
    Ok(())
}

#[test]
fn ensure_plain_bed_windows_not_empty_errors_when_no_windows_survive() -> Result<()> {
    // Arrange
    let mut windows_map = FxHashMap::default();
    windows_map.insert("chr1".to_string(), Windows::from_tuples(&[])?);
    windows_map.insert("chr2".to_string(), Windows::from_tuples(&[])?);

    // Act
    let result = ensure_plain_bed_windows_not_empty(&windows_map);

    // Assert
    let err = result.expect_err("empty selected BED windows should error");
    assert!(
        err.to_string()
            .contains("BED file did not contain any valid windows on the selected chromosomes")
    );
    Ok(())
}

#[test]
fn ensure_plain_bed_windows_not_empty_accepts_one_surviving_window() -> Result<()> {
    // Arrange
    let mut windows_map = FxHashMap::default();
    windows_map.insert("chr1".to_string(), Windows::from_tuples(&[])?);
    windows_map.insert("chr2".to_string(), Windows::from_tuples(&[(10, 20, 0)])?);

    // Act / Assert
    ensure_plain_bed_windows_not_empty(&windows_map)?;
    Ok(())
}

#[test]
fn ensure_grouped_bed_windows_not_empty_errors_when_no_grouped_windows_survive() -> Result<()> {
    // Arrange
    let mut windows_map = FxHashMap::default();
    windows_map.insert("chr1".to_string(), GroupedWindows::from_tuples(&[], None)?);
    windows_map.insert("chr2".to_string(), GroupedWindows::from_tuples(&[], None)?);

    // Act
    let result = ensure_grouped_bed_windows_not_empty(&windows_map);

    // Assert
    let err = result.expect_err("empty selected grouped BED windows should error");
    assert!(err.to_string().contains(
        "grouped BED file did not contain any valid windows on the selected chromosomes"
    ));
    Ok(())
}

#[test]
fn ensure_grouped_bed_windows_not_empty_accepts_one_surviving_grouped_window() -> Result<()> {
    // Arrange
    let mut windows_map = FxHashMap::default();
    windows_map.insert("chr1".to_string(), GroupedWindows::from_tuples(&[], None)?);
    windows_map.insert(
        "chr2".to_string(),
        GroupedWindows::from_tuples(&[(10, 20, 0)], None)?,
    );

    // Act / Assert
    ensure_grouped_bed_windows_not_empty(&windows_map)?;
    Ok(())
}

#[test]
fn build_bin_info_uses_size_offsets_in_output_indices() -> Result<()> {
    // Arrange
    let chromosomes = vec!["chr1".to_string(), "chr2".to_string()];
    let contigs = contigs_with_lengths(&[("chr1", 150), ("chr2", 80)]);
    let blacklist_map = FxHashMap::default();
    let mut chr_offsets = FxHashMap::default();
    chr_offsets.insert("chr1".to_string(), 0);
    chr_offsets.insert("chr2".to_string(), 2);

    // Act
    let bins = build_bin_info(
        &WindowSpec::Size(100),
        &chromosomes,
        &contigs,
        None,
        &blacklist_map,
        &chr_offsets,
    )?;

    // Assert
    assert_eq!(
        bins[0],
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 0,
            end: 100,
            output_index: 0,
            blacklisted_fraction: 0.0,
        }
    );
    assert_eq!(
        bins[1],
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 100,
            end: 150,
            output_index: 1,
            blacklisted_fraction: 0.0,
        }
    );
    assert_eq!(
        bins[2],
        WindowBinInfo {
            chromosome: "chr2".to_string(),
            start: 0,
            end: 80,
            output_index: 2,
            blacklisted_fraction: 0.0,
        }
    );
    Ok(())
}

#[test]
fn build_bin_info_preserves_bed_original_indices_and_sorts_by_them() -> Result<()> {
    // Arrange
    let chromosomes = vec!["chr1".to_string()];
    let contigs = contigs_with_lengths(&[("chr1", 100)]);
    let mut windows_map = FxHashMap::default();
    windows_map.insert(
        "chr1".to_string(),
        Windows::from_tuples(&[(20, 30, 9_u64), (0, 10, 3_u64)])?,
    );
    let blacklist_map = FxHashMap::default();
    let chr_offsets = FxHashMap::default();

    // Act
    let bins = build_bin_info(
        &WindowSpec::Bed(PathBuf::from("windows.bed")),
        &chromosomes,
        &contigs,
        Some(&windows_map),
        &blacklist_map,
        &chr_offsets,
    )?;

    // Assert
    assert_eq!(
        bins[0],
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 0,
            end: 10,
            output_index: 3,
            blacklisted_fraction: 0.0,
        }
    );
    assert_eq!(
        bins[1],
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 20,
            end: 30,
            output_index: 9,
            blacklisted_fraction: 0.0,
        }
    );
    Ok(())
}
