use super::*;

fn grouped_windows(entries: &[(u64, u64, u64)]) -> GroupedWindows {
    GroupedWindows::from_tuples(entries, None).expect("test grouped windows should be valid")
}

fn group_names(entries: &[(u64, &str)]) -> FxHashMap<u64, String> {
    let mut names = FxHashMap::default();
    for (group_idx, group_name) in entries {
        names.insert(*group_idx, (*group_name).to_string());
    }
    names
}

fn layout_segments_for_chr(
    layout: &GroupedCoverageLayout,
    chromosome: &str,
) -> Vec<(u64, u64, u64)> {
    layout
        .segments_by_chr
        .get(chromosome)
        .map(|windows| {
            windows
                .as_slice()
                .iter()
                .map(|segment| (segment.start(), segment.end(), segment.idx()))
                .collect()
        })
        .unwrap_or_default()
}

#[test]
fn should_error_when_selected_column_6_later_contains_invalid_strand() -> Result<()> {
    // Arrange:
    // - The sampled data rows make column 6 the detected strand column.
    // - The next row then contains an invalid strand token in that already-selected column.
    // - Detection is intentionally bounded, but parsing must stay strict for the full file.
    let mut bed = tempfile::NamedTempFile::new()?;
    for row_idx in 0..GROUPED_BED_STRAND_SAMPLE_ROWS {
        writeln!(
            bed,
            "chr1\t{}\t{}\talpha\t0\t+",
            row_idx * 10,
            row_idx * 10 + 5
        )?;
    }
    writeln!(
        bed,
        "chr1\t{}\t{}\talpha\t0\tx",
        GROUPED_BED_STRAND_SAMPLE_ROWS * 10,
        GROUPED_BED_STRAND_SAMPLE_ROWS * 10 + 5
    )?;

    // Act
    let error = load_grouped_windows_from_bed(
        bed.path(),
        None,
        true,
        None,
        Some((GROUPED_BED_STRAND_SAMPLE_ROWS + 1) as u64),
    )
    .expect_err("invalid strand after the sampling window should fail during parsing");

    // Assert
    assert!(
        error.to_string().contains("invalid strand 'x' in column 6"),
        "unexpected error: {error:?}"
    );
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_keeps_raw_segments_when_unique_bases_is_disabled() -> Result<()> {
    // Arrange
    // Group 0 has two overlapping windows:
    //   [10, 20) and [15, 25)
    // In raw mode, both must remain as separate segments and their lengths must both
    // contribute to `group_span_positions`
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert(
        "chr1".to_string(),
        grouped_windows(&[(10, 20, 0), (15, 25, 0), (40, 45, 1)]),
    );
    let group_idx_to_name = group_names(&[(0, "alpha"), (1, "beta")]);

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_idx_to_name,
        &["chr1".to_string()],
        false,
    )?;

    // Assert
    assert_eq!(
        layout_segments_for_chr(&layout, "chr1"),
        vec![(10, 20, 0), (15, 25, 1), (40, 45, 2)]
    );
    assert_eq!(layout.segment_idx_to_group_idx.get(&0), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&1), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&2), Some(&1));
    assert_eq!(layout.group_span_positions.get(&0), Some(&20));
    assert_eq!(layout.group_span_positions.get(&1), Some(&5));
    assert_eq!(
        layout.group_idx_to_name.get(&0).map(String::as_str),
        Some("alpha")
    );
    assert_eq!(
        layout.group_idx_to_name.get(&1).map(String::as_str),
        Some("beta")
    );
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_merges_same_group_overlaps_touches_and_duplicates() -> Result<()> {
    // Arrange
    // Group 0 contributes:
    //   [10, 20) overlaps [15, 25)
    //   [25, 30) touches the merged tail at 25
    // So all three must collapse to [10, 30)
    //
    // Group 1 contributes:
    //   [40, 45) and an identical duplicate [40, 45)
    // Duplicates also collapse to one merged segment
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert(
        "chr1".to_string(),
        grouped_windows(&[
            (10, 20, 0),
            (15, 25, 0),
            (25, 30, 0),
            (40, 45, 1),
            (40, 45, 1),
        ]),
    );
    let group_idx_to_name = group_names(&[(0, "alpha"), (1, "beta")]);

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_idx_to_name,
        &["chr1".to_string()],
        true,
    )?;

    // Assert
    assert_eq!(
        layout_segments_for_chr(&layout, "chr1"),
        vec![(10, 30, 0), (40, 45, 1)]
    );
    assert_eq!(layout.segment_idx_to_group_idx.get(&0), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&1), Some(&1));
    assert_eq!(layout.group_span_positions.get(&0), Some(&20));
    assert_eq!(layout.group_span_positions.get(&1), Some(&5));
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_does_not_merge_overlaps_from_different_groups() -> Result<()> {
    // Arrange
    // Group 0 covers [10, 20)
    // Group 1 covers [15, 25) and [25, 30), which merge within group 1 to [15, 30)
    // The cross-group overlap [15, 20) must remain represented in both groups
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert(
        "chr1".to_string(),
        grouped_windows(&[(10, 20, 0), (15, 25, 1), (25, 30, 1)]),
    );

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_names(&[(0, "alpha"), (1, "beta")]),
        &["chr1".to_string()],
        true,
    )?;

    // Assert
    assert_eq!(
        layout_segments_for_chr(&layout, "chr1"),
        vec![(10, 20, 0), (15, 30, 1)]
    );
    assert_eq!(layout.segment_idx_to_group_idx.get(&0), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&1), Some(&1));
    assert_eq!(layout.group_span_positions.get(&0), Some(&10));
    assert_eq!(layout.group_span_positions.get(&1), Some(&15));
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_keeps_disjoint_segments_for_the_same_group() -> Result<()> {
    // Arrange
    // Group 0 has two separated spans, so unique-base mode must keep two segments
    // rather than inventing one larger interval with a gap in the middle
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert(
        "chr1".to_string(),
        grouped_windows(&[(0, 5, 0), (10, 15, 0)]),
    );

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_names(&[(0, "alpha")]),
        &["chr1".to_string()],
        true,
    )?;

    // Assert
    assert_eq!(
        layout_segments_for_chr(&layout, "chr1"),
        vec![(0, 5, 0), (10, 15, 1)]
    );
    assert_eq!(layout.segment_idx_to_group_idx.get(&0), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&1), Some(&0));
    assert_eq!(layout.group_span_positions.get(&0), Some(&10));
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_assigns_segment_indices_by_requested_chromosome_order()
-> Result<()> {
    // Arrange
    // The requested chromosome order is chr3, chr1, chr2.
    // chr2 is missing and must be skipped without consuming an index.
    // That means chr3 gets segment_idx 0 and chr1 gets segment_idx 1
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert("chr1".to_string(), grouped_windows(&[(20, 30, 1)]));
    grouped_windows_by_chr.insert("chr3".to_string(), grouped_windows(&[(5, 10, 0)]));

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_names(&[(0, "alpha"), (1, "beta")]),
        &["chr3".to_string(), "chr2".to_string(), "chr1".to_string()],
        false,
    )?;

    // Assert
    assert_eq!(layout_segments_for_chr(&layout, "chr3"), vec![(5, 10, 0)]);
    assert_eq!(layout_segments_for_chr(&layout, "chr1"), vec![(20, 30, 1)]);
    assert!(layout_segments_for_chr(&layout, "chr2").is_empty());
    assert_eq!(layout.segment_idx_to_group_idx.get(&0), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&1), Some(&1));
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_sorts_segments_by_start_end_and_group_idx_within_chromosome()
-> Result<()> {
    // Arrange
    // After merging in unique-base mode, the chromosome contains three segments with the same
    // start coordinate:
    //   group 1 -> [10, 15)
    //   group 0 -> [10, 20)
    //   group 2 -> [12, 18)
    //
    // Sorting by (start, end, group_idx) must therefore yield:
    //   [10, 15) group 1
    //   [10, 20) group 0
    //   [12, 18) group 2
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert(
        "chr1".to_string(),
        grouped_windows(&[(10, 20, 0), (10, 15, 1), (12, 18, 2)]),
    );

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_names(&[(0, "alpha"), (1, "beta"), (2, "gamma")]),
        &["chr1".to_string()],
        true,
    )?;

    // Assert
    assert_eq!(
        layout_segments_for_chr(&layout, "chr1"),
        vec![(10, 15, 0), (10, 20, 1), (12, 18, 2)]
    );
    assert_eq!(layout.segment_idx_to_group_idx.get(&0), Some(&1));
    assert_eq!(layout.segment_idx_to_group_idx.get(&1), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&2), Some(&2));
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_returns_empty_segments_but_preserves_group_names_for_empty_input()
-> Result<()> {
    // Arrange
    // Even with no grouped windows, the caller-provided group-name map should survive unchanged
    // so downstream metadata can still be written deterministically
    let grouped_windows_by_chr = FxHashMap::default();
    let group_idx_to_name = group_names(&[(0, "alpha"), (1, "beta")]);

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_idx_to_name,
        &["chr1".to_string()],
        true,
    )?;

    // Assert
    assert!(layout.segments_by_chr.is_empty());
    assert!(layout.segment_idx_to_group_idx.is_empty());
    assert!(layout.group_span_positions.is_empty());
    assert_eq!(layout.group_idx_to_name, group_idx_to_name);
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_preserves_unused_group_names_in_nonempty_layout() -> Result<()> {
    // Arrange
    // Group 2 is present in the metadata map but has no windows in the layout.
    // The layout builder should not delete that name, because the caller may still need the
    // full sidecar mapping that came from the grouped BED loader
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert("chr1".to_string(), grouped_windows(&[(0, 10, 0)]));
    let group_idx_to_name = group_names(&[(0, "alpha"), (2, "gamma")]);

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_idx_to_name,
        &["chr1".to_string()],
        false,
    )?;

    // Assert
    assert_eq!(layout_segments_for_chr(&layout, "chr1"), vec![(0, 10, 0)]);
    assert_eq!(layout.group_idx_to_name, group_idx_to_name);
    assert_eq!(layout.group_span_positions.get(&0), Some(&10));
    assert!(layout.group_span_positions.get(&2).is_none());
    Ok(())
}

#[test]
fn build_grouped_coverage_layout_keeps_chromosomes_independent_while_summing_group_spans()
-> Result<()> {
    // Arrange
    // Group 0 appears on two chromosomes:
    // - chr1: [0, 5) and [5, 10) touch, so unique-base mode merges them to [0, 10)
    // - chr2: [20, 25) stays separate because chromosomes must never merge into each other
    //
    // The per-group span should therefore be 10 + 5 = 15 across the two chromosomes
    let mut grouped_windows_by_chr = FxHashMap::default();
    grouped_windows_by_chr.insert(
        "chr1".to_string(),
        grouped_windows(&[(0, 5, 0), (5, 10, 0)]),
    );
    grouped_windows_by_chr.insert("chr2".to_string(), grouped_windows(&[(20, 25, 0)]));

    // Act
    let layout = build_grouped_coverage_layout(
        &grouped_windows_by_chr,
        &group_names(&[(0, "alpha")]),
        &["chr1".to_string(), "chr2".to_string()],
        true,
    )?;

    // Assert
    assert_eq!(layout_segments_for_chr(&layout, "chr1"), vec![(0, 10, 0)]);
    assert_eq!(layout_segments_for_chr(&layout, "chr2"), vec![(20, 25, 1)]);
    assert_eq!(layout.segment_idx_to_group_idx.get(&0), Some(&0));
    assert_eq!(layout.segment_idx_to_group_idx.get(&1), Some(&0));
    assert_eq!(layout.group_span_positions.get(&0), Some(&15));
    Ok(())
}
