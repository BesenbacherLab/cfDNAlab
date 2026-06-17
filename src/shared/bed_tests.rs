use super::*;

mod tests_newest_bed_loaders {

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
    fn build_grouped_coverage_layout_keeps_raw_segments_when_unique_bases_is_disabled() -> Result<()>
    {
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
    fn build_grouped_coverage_layout_merges_same_group_overlaps_touches_and_duplicates()
    -> Result<()> {
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
    fn build_grouped_coverage_layout_preserves_unused_group_names_in_nonempty_layout() -> Result<()>
    {
        // Arrange
        // Group 2 is present in the metadata map but has no windows in the layout.
        // The layout builder should not delete that name, because the caller may still need the
        // full metadata mapping that came from the grouped BED loader
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
}

mod tests_bed_loader {
    use crate::shared::bed::{
        GroupedBedStrandColumn, Strand, load_grouped_windows_from_bed,
        load_scored_windows_from_bed, load_windows_from_bed, write_group_idx_to_name_tsv,
    };
    use crate::shared::interval::{IndexedInterval, ScoredInterval};
    use anyhow::Result;
    use flate2::{Compression, write::GzEncoder};
    use fxhash::FxHashMap;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    /// Write helper content to a temporary BED file used across tests.
    fn write_bed(lines: &[&str]) -> Result<NamedTempFile> {
        let mut file = NamedTempFile::new()?;
        for line in lines {
            writeln!(file, "{}", line)?;
        }
        Ok(file)
    }

    fn indexed_windows(entries: &[(u64, u64, u64)]) -> Vec<IndexedInterval<u64>> {
        entries
            .iter()
            .map(|&(start, end, original_index)| {
                IndexedInterval::new(start, end, original_index)
                    .expect("test windows should be valid non-empty intervals")
            })
            .collect()
    }

    fn scored_windows(entries: &[(u64, u64, u64, f64)]) -> Vec<ScoredInterval<u64>> {
        entries
            .iter()
            .map(|&(start, end, original_index, score)| {
                ScoredInterval::new(start, end, original_index, score)
                    .expect("test scored windows should be valid non-empty intervals")
            })
            .collect()
    }

    #[test]
    fn should_keep_only_whitelisted_chromosomes_when_loading_bed() -> Result<()> {
        // Arrange
        let bed = write_bed(&["chr1\t0\t10", "chr2\t5\t15", "chr1\t20\t30"])?;
        let whitelist = vec!["chr1".to_string()];

        // Act
        let map = load_windows_from_bed(bed.path(), Some(whitelist.as_slice()), None, None)?;

        // Assert
        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(chr1.as_slice(), indexed_windows(&[(0, 10, 0), (20, 30, 2)]));

        let empty = load_windows_from_bed(
            bed.path(),
            Some(["chr3".to_string()].as_slice()),
            None,
            None,
        )?;
        assert!(
            empty
                .get("chr3")
                .expect("chr3 entry should exist")
                .as_slice()
                .is_empty()
        );
        Ok(())
    }

    #[test]
    fn should_filter_windows_by_predicate_when_loading_bed() -> Result<()> {
        // Arrange
        let bed = write_bed(&["chr1\t0\t5", "chr1\t10\t25", "chr1\t30\t33"])?;
        let keep_large = |_: &str, start: u64, end: u64| (end - start) >= 10;

        // Act
        let map = load_windows_from_bed(bed.path(), None, Some(&keep_large), None)?;

        // Assert
        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(chr1.as_slice(), indexed_windows(&[(10, 25, 1)]));
        Ok(())
    }

    #[test]
    fn should_load_gzipped_bed() -> Result<()> {
        let gz = tempfile::Builder::new().suffix(".bed.gz").tempfile()?;

        {
            let file = std::fs::File::create(gz.path())?;
            let mut encoder = GzEncoder::new(file, Compression::default());
            writeln!(encoder, "chr1\t0\t5")?;
            writeln!(encoder, "chr1\t10\t15")?;
            encoder.finish()?;
        }

        let map = load_windows_from_bed(gz.path(), None, None, None)?;
        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(chr1.as_slice(), indexed_windows(&[(0, 5, 0), (10, 15, 1)]));
        Ok(())
    }

    #[test]
    fn should_validate_expected_window_count_with_whitelist() -> Result<()> {
        // Arrange
        let bed = write_bed(&["chr1\t0\t4", "chr2\t4\t8", "chr2\t8\t12"])?;
        let whitelist = vec!["chr2".to_string()];

        // Act
        let map = load_windows_from_bed(bed.path(), Some(whitelist.as_slice()), None, Some(3))?;

        // Assert: only the allowed chromosome is returned, but **original indices include skipped windows**
        let chr2 = map.get("chr2").expect("chr2 entry missing");
        assert_eq!(chr2.as_slice(), indexed_windows(&[(4, 8, 1), (8, 12, 2)]));

        // And mismatched expectations yield an error
        let err = load_windows_from_bed(bed.path(), Some(whitelist.as_slice()), None, Some(2))
            .expect_err("expected incorrect exp_num_windows to error");
        assert!(
            err.to_string()
                .contains("did not contain the correct number of windows"),
            "unexpected error: {err:?}"
        );
        Ok(())
    }

    #[test]
    fn should_error_on_invalid_windows_even_with_expected_count() -> Result<()> {
        // Arrange: second line has end <= start, so it should error regardless of exp_num_windows
        let bed = write_bed(&["chr1\t0\t5", "chr1\t5\t4", "chr2\t10\t20"])?;

        // Act + Assert
        let err = load_windows_from_bed(bed.path(), None, None, Some(3))
            .expect_err("invalid window should fail loading");
        assert!(
            err.to_string()
                .contains("end (4) must be greater than start (5)"),
            "unexpected error: {err:?}"
        );
        Ok(())
    }

    #[test]
    fn should_sort_grouped_windows_and_reuse_group_indices_when_loading_bed() -> Result<()> {
        // Arrange:
        // - Group indices are assigned when each group name is first seen in file order.
        // - "beta" appears first, so it must get group_idx 0.
        // - "alpha" appears second, so it must get group_idx 1.
        // - Within each chromosome, grouped windows are then sorted by (start, end), but that sort
        //   must not rewrite the already-assigned group indices.
        let bed = write_bed(&[
            "chr2\t20\t30\tbeta",
            "chr1\t15\t18\tbeta",
            "chr1\t10\t12\talpha",
            "chr2\t5\t8\talpha",
        ])?;

        let (map, group_idx_to_name, strand_detection) =
            load_grouped_windows_from_bed(bed.path(), None, false, None, Some(4))?;
        assert!(
            strand_detection.is_none(),
            "strand detection should not run when read_strands is false"
        );

        assert_eq!(map.len(), 2);
        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(
            chr1.windows_as_slice(),
            indexed_windows(&[(10, 12, 1), (15, 18, 0)])
        );
        assert!(
            chr1.strands.is_none(),
            "strands should not be loaded when read_strands is false"
        );

        let chr2 = map.get("chr2").expect("chr2 missing");
        assert_eq!(
            chr2.windows_as_slice(),
            indexed_windows(&[(5, 8, 1), (20, 30, 0)])
        );

        assert_eq!(group_idx_to_name.len(), 2);
        assert_eq!(group_idx_to_name.get(&0).map(String::as_str), Some("beta"));
        assert_eq!(group_idx_to_name.get(&1).map(String::as_str), Some("alpha"));
        Ok(())
    }

    #[test]
    fn should_read_grouped_bed_strands_from_column_6() -> Result<()> {
        // Arrange:
        // - Column 4 is the group name, column 5 is a BED score-like value, and column 6 is strand.
        // - The loader sorts by coordinate, so strand values must move with their source rows.
        let bed = write_bed(&[
            "chr1\t20\t30\tbeta\t0\t-",
            "chr1\t10\t15\talpha\t0\t+",
            "chr1\t15\t18\tbeta\t0\t.",
        ])?;

        let (map, _group_idx_to_name, strand_detection) =
            load_grouped_windows_from_bed(bed.path(), None, true, None, Some(3))?;
        let strand_detection = strand_detection.expect("strand detection should run");
        assert_eq!(
            strand_detection.column,
            Some(GroupedBedStrandColumn::Column6)
        );

        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(
            chr1.windows_as_slice(),
            indexed_windows(&[(10, 15, 1), (15, 18, 0), (20, 30, 0)])
        );
        assert_eq!(
            chr1.strands
                .as_ref()
                .expect("strands should be loaded")
                .as_slice(),
            &[Strand::Forward, Strand::Unstranded, Strand::Reverse]
        );
        Ok(())
    }

    #[test]
    fn should_read_grouped_bed_strands_from_column_5_when_no_column_6_exists() -> Result<()> {
        // Arrange:
        // - This is a non-standard grouped BED shape where column 4 is the group and column 5 is
        //   strand. It is only accepted because there is no column 6.
        let bed = write_bed(&["chr1\t10\t15\talpha\t+", "chr1\t20\t25\tbeta\t-"])?;

        let (map, _group_idx_to_name, strand_detection) =
            load_grouped_windows_from_bed(bed.path(), None, true, None, Some(2))?;
        let strand_detection = strand_detection.expect("strand detection should run");
        assert_eq!(
            strand_detection.column,
            Some(GroupedBedStrandColumn::Column5)
        );

        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(
            chr1.windows_as_slice(),
            indexed_windows(&[(10, 15, 0), (20, 25, 1)])
        );
        assert_eq!(
            chr1.strands
                .as_ref()
                .expect("strands should be loaded")
                .as_slice(),
            &[Strand::Forward, Strand::Reverse]
        );
        Ok(())
    }

    #[test]
    fn should_error_when_column_5_looks_stranded_but_column_6_exists_without_strands() -> Result<()>
    {
        // Arrange:
        // - With 6 columns, strand belongs in column 6.
        // - A strand-looking column 5 is ambiguous because it could be a non-standard file or a
        //   misplaced strand column, so the loader must not silently treat the file as unstranded.
        let bed = write_bed(&["chr1\t10\t15\talpha\t+\t0", "chr1\t20\t25\tbeta\t-\t0"])?;

        let error = load_grouped_windows_from_bed(bed.path(), None, true, None, Some(2))
            .expect_err("ambiguous strand columns should fail");

        assert!(
            error
                .to_string()
                .contains("When 6 or more BED columns are supplied, put strands in column 6"),
            "unexpected error: {error:?}"
        );
        Ok(())
    }

    #[test]
    fn should_treat_wide_grouped_bed_as_unstranded_when_no_strand_column_is_detected() -> Result<()>
    {
        // Arrange:
        // - The file has 6 columns, but neither column 5 nor column 6 contains UCSC strand tokens.
        // - This is a wide non-standard grouped BED-like file, so the loader should keep the intervals
        //   and report that no strand column was selected.
        let bed = write_bed(&[
            "chr1\t20\t25\tbeta\t7.2\tannotation_b",
            "chr1\t10\t15\talpha\t3.1\tannotation_a",
        ])?;

        // Act
        let (map, _group_idx_to_name, strand_detection) =
            load_grouped_windows_from_bed(bed.path(), None, true, None, Some(2))?;

        // Assert
        let strand_detection = strand_detection.expect("strand detection should run");
        assert_eq!(strand_detection.column, None);
        assert!(
            strand_detection.saw_column6,
            "detection metadata should record that a wide file was sampled"
        );

        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(
            chr1.windows_as_slice(),
            indexed_windows(&[(10, 15, 1), (20, 25, 0)])
        );
        assert!(
            chr1.strands.is_none(),
            "wide files without strand tokens in the sampled strand columns should be unstranded"
        );
        Ok(())
    }

    #[test]
    fn should_filter_scored_windows_but_preserve_original_indices_and_sorting() -> Result<()> {
        // Arrange:
        // - Original line indices are 0:[30,40) score 0.5, 1:[10,15) score 2.0, 2:[20,25) score 3.5.
        // - Filtering keeps only scores >= 2.0, so lines 1 and 2 survive.
        // - The retained output must then be sorted by genomic coordinates, giving [10,15) before
        //   [20,25), while preserving their original indices 1 and 2.
        // - Span should therefore be min start 10 and max end 25.
        let bed = write_bed(&[
            "chr1\t30\t40\t0.5",
            "chr1\t10\t15\t2.0",
            "chr1\t20\t25\t3.5",
        ])?;
        let keep_high_scores = |_: &str, _: u64, _: u64, score: f64| score >= 2.0;

        let map = load_scored_windows_from_bed(bed.path(), None, Some(&keep_high_scores), Some(3))?;

        assert_eq!(map.len(), 1);
        let chr1 = map.get("chr1").expect("chr1 missing");
        assert_eq!(
            chr1.as_slice(),
            scored_windows(&[(10, 15, 1, 2.0), (20, 25, 2, 3.5)])
        );
        assert_eq!(chr1.span_start(), 10);
        assert_eq!(chr1.span_end(), 25);
        Ok(())
    }

    #[test]
    fn should_keep_only_whitelisted_chromosomes_when_loading_grouped_bed() -> Result<()> {
        // Arrange:
        // - Whitelist keeps only chr2, but group indices are assigned in original file order.
        // - "alpha" appears first on chr1 -> group_idx 0 would be tempting, but the grouped loader
        //   only assigns indices when a kept row is processed.
        // - The kept chr2 rows appear as "beta" first and then "alpha", so chr2 must contain
        //   [10,15) -> 0 and [5,9) -> 1, later sorted into genomic order.
        let bed = write_bed(&[
            "chr1\t0\t5\talpha",
            "chr2\t10\t15\tbeta",
            "chr2\t5\t9\talpha",
        ])?;
        let whitelist = vec!["chr2".to_string()];

        let (map, group_idx_to_name, _strand_detection) = load_grouped_windows_from_bed(
            bed.path(),
            Some(whitelist.as_slice()),
            false,
            None,
            Some(3),
        )?;

        assert_eq!(map.len(), 1);
        assert!(
            map.get("chr1").is_none(),
            "chr1 should be excluded by the chromosome whitelist"
        );
        let chr2 = map.get("chr2").expect("chr2 missing");
        assert_eq!(
            chr2.windows_as_slice(),
            indexed_windows(&[(5, 9, 1), (10, 15, 0)])
        );
        assert_eq!(group_idx_to_name.get(&0).map(String::as_str), Some("beta"));
        assert_eq!(group_idx_to_name.get(&1).map(String::as_str), Some("alpha"));
        Ok(())
    }

    #[test]
    fn should_error_when_grouped_bed_is_missing_group_name() -> Result<()> {
        let bed = write_bed(&["chr1\t0\t10"])?;

        let error = load_grouped_windows_from_bed(bed.path(), None, false, None, None)
            .expect_err("missing group name should fail");

        assert!(error.to_string().contains("missing group name"));
        Ok(())
    }

    #[test]
    fn should_keep_original_indices_when_loading_scored_bed_with_whitelist() -> Result<()> {
        // Arrange:
        // - File order gives original indices 0:[0,4), 1:[9,12), 2:[4,8).
        // - Whitelisting chr2 removes line 0 but must keep the surviving original indices 1 and 2.
        // - Sorting within chr2 then places [4,8) before [9,12), so the expected order is idx 2 then 1.
        let bed = write_bed(&["chr1\t0\t4\t1.0", "chr2\t9\t12\t2.0", "chr2\t4\t8\t3.0"])?;
        let whitelist = vec!["chr2".to_string()];

        let map =
            load_scored_windows_from_bed(bed.path(), Some(whitelist.as_slice()), None, Some(3))?;

        assert_eq!(map.len(), 1);
        assert!(
            map.get("chr1").is_none(),
            "chr1 should be excluded by the chromosome whitelist"
        );
        let chr2 = map.get("chr2").expect("chr2 missing");
        assert_eq!(
            chr2.as_slice(),
            scored_windows(&[(4, 8, 2, 3.0), (9, 12, 1, 2.0)])
        );
        assert_eq!(chr2.span_start(), 4);
        assert_eq!(chr2.span_end(), 12);
        Ok(())
    }

    #[test]
    fn should_error_when_scored_bed_is_missing_score() -> Result<()> {
        let bed = write_bed(&["chr1\t0\t10"])?;

        let error = load_scored_windows_from_bed(bed.path(), None, None, None)
            .expect_err("missing score should fail");

        assert!(error.to_string().contains("missing score"));
        Ok(())
    }

    #[test]
    fn should_error_when_scored_bed_has_invalid_score() -> Result<()> {
        let bed = write_bed(&["chr1\t0\t10\tnot_a_float"])?;

        let error = load_scored_windows_from_bed(bed.path(), None, None, None)
            .expect_err("invalid score should fail");

        assert!(error.to_string().contains("invalid score 'not_a_float'"));
        Ok(())
    }

    #[test]
    fn should_write_group_index_tsv_sorted_and_sanitized() -> Result<()> {
        // Arrange:
        // - Output rows are written in increasing numeric group index order, so 0 must precede 2.
        // - Embedded newlines are replaced with spaces, and tabs are expanded so the TSV stays
        //   one logical row per group.
        let temp = TempDir::new()?;
        let path = temp.path().join("group_index.tsv");
        let mut group_idx_to_name = FxHashMap::default();
        group_idx_to_name.insert(2_u64, "beta\tname".to_string());
        group_idx_to_name.insert(0_u64, "alpha\nname".to_string());

        write_group_idx_to_name_tsv(&path, &group_idx_to_name)?;

        let written = std::fs::read_to_string(&path)?;
        let lines: Vec<_> = written.lines().collect();
        assert_eq!(lines[0], "group_idx\tgroup_name");
        assert_eq!(lines[1], "0\talpha name");
        assert_eq!(lines[2], "2\tbeta    name");
        Ok(())
    }

    #[cfg(feature = "cmd_prepare_windows")]
    #[test]
    fn should_detect_header_after_comments_and_blank_lines() -> Result<()> {
        use crate::shared::bed::detect_header;

        // The detector should skip comments and empty lines, so the first meaningful line here is the
        // literal header "chrom\tstart\tend".
        let bed = write_bed(&["# comment", "", "chrom\tstart\tend", "chr1\t0\t10"])?;

        assert!(detect_header(bed.path(), '\t')?);
        Ok(())
    }

    #[cfg(feature = "cmd_prepare_windows")]
    #[test]
    fn should_detect_coordinate_lines_without_header() {
        use crate::shared::bed::line_looks_like_header;

        // "chr1\t0\t10" has numeric coordinate columns and is therefore data, not a header.
        // The literal column names and comment lines should be treated as header-like.
        assert!(!line_looks_like_header("chr1\t0\t10", '\t'));
        assert!(line_looks_like_header("chrom\tstart\tend", '\t'));
        assert!(line_looks_like_header("# comment", '\t'));
    }
}

mod tests_flattening {
    use crate::shared::bed::*;
    use crate::shared::interval::{IndexedInterval, ScoredInterval, Span};

    // Helper: build a start-sorted Windows from (s,e) pairs (original_idx is dummy)
    fn mk_sorted(pairs: &[(u64, u64)]) -> Windows {
        let windows = pairs
            .iter()
            .enumerate()
            .map(|(i, &(start, end))| {
                IndexedInterval::new(start, end, i as u64).expect("test windows should be valid")
            })
            .collect();
        Windows::from_sorted(windows)
    }

    // Helper: assert strictly sorted and non-overlapping (touching should have been merged away)
    fn assert_sorted_non_overlapping(ws: &[IndexedInterval<u64>]) {
        for index in 1..ws.len() {
            assert!(
                ws[index - 1].start() <= ws[index].start(),
                "not sorted: prev.start={} > cur.start={}",
                ws[index - 1].start(),
                ws[index].start()
            );
            assert!(
                ws[index - 1].end() < ws[index].start(),
                "intervals overlap or still touch: prev={:?}, cur={:?}",
                ws[index - 1],
                ws[index]
            );
            assert!(
                ws[index - 1].start() < ws[index - 1].end(),
                "invalid interval with zero/negative length"
            );
        }
        if let Some(last) = ws.last() {
            assert!(
                last.start() < last.end(),
                "invalid interval with zero/negative length"
            );
        }
    }

    // Helper: assert indices are sequential starting at start_idx
    fn assert_sequential_indices(ws: &[IndexedInterval<u64>], start_idx: u64) {
        for (index, window) in ws.iter().enumerate() {
            assert_eq!(
                window.idx(),
                start_idx + index as u64,
                "non-sequential index at k={}",
                index
            );
        }
    }

    #[test]
    fn flatten_empty() {
        let w = Windows::from_sorted(Vec::new());
        let (flat, next) = w.into_flattened_reindexed(0);
        assert!(flat.as_slice().is_empty());
        assert_eq!(next, 0);
        assert_eq!(flat.span_start(), 0);
        assert_eq!(flat.span_end(), 0);
    }

    #[test]
    fn flatten_singleton() {
        let w = mk_sorted(&[(10, 20)]);
        let (flat, next) = w.into_flattened_reindexed(7);
        let a = flat.as_slice();
        assert_eq!(a[0].into_tuple(), (10, 20, 7));
        assert_eq!(next, 8);
        assert_eq!(flat.span_start(), 10);
        assert_eq!(flat.span_end(), 20);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 7);
    }

    #[test]
    fn flatten_no_merges() {
        // Non-overlapping and non-touching
        let w = mk_sorted(&[(10, 15), (20, 25), (30, 35)]);
        let (flat, next) = w.into_flattened_reindexed(0);
        let a = flat.as_slice();
        assert_eq!(a.len(), 3);
        assert_eq!(next, 3);
        // Starts/ends preserved
        assert_eq!(a[0].start(), 10);
        assert_eq!(a[0].end(), 15);
        assert_eq!(a[1].start(), 20);
        assert_eq!(a[1].end(), 25);
        assert_eq!(a[2].start(), 30);
        assert_eq!(a[2].end(), 35);
        assert_eq!(flat.span_start(), 10);
        assert_eq!(flat.span_end(), 35);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 0);
    }

    #[test]
    fn flatten_touching_merges() {
        // Touching intervals must merge (half-open semantics)
        let w = mk_sorted(&[(10, 15), (15, 20), (20, 30)]);
        let (flat, next) = w.into_flattened_reindexed(100);
        let a = flat.as_slice();
        assert_eq!(a[0].into_tuple(), (10, 30, 100));
        assert_eq!(next, 101);
        assert_eq!(flat.span_start(), 10);
        assert_eq!(flat.span_end(), 30);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 100);
    }

    #[test]
    fn flatten_overlapping_chain() {
        // Mixed: one disjoint small block and a chain that overlaps/touches
        let w = mk_sorted(&[(5, 7), (10, 14), (12, 16), (16, 19)]);
        let (flat, next) = w.into_flattened_reindexed(50);
        let a = flat.as_slice();
        // Expect (5,7) and (10,19)
        assert_eq!(a.len(), 2);
        assert_eq!(a[0].into_tuple(), (5, 7, 50));
        assert_eq!(a[1].into_tuple(), (10, 19, 51));
        assert_eq!(next, 52);
        assert_eq!(flat.span_start(), 5);
        assert_eq!(flat.span_end(), 19);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 50);
    }

    #[test]
    fn flatten_large_start_idx() {
        // Sanity: indices carry forward correctly from large start
        let w = mk_sorted(&[(0, 1), (2, 3), (4, 5)]);
        let (flat, next) = w.into_flattened_reindexed(1_000_000);
        let a = flat.as_slice();
        assert_eq!(a.len(), 3);
        assert_eq!(a[0].idx(), 1_000_000);
        assert_eq!(a[1].idx(), 1_000_001);
        assert_eq!(a[2].idx(), 1_000_002);
        assert_eq!(next, 1_000_003);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 1_000_000);
    }

    #[test]
    fn grouped_windows_sort_and_preserve_group_indices() {
        // Arrange:
        // - Inputs are unsorted by start, but group indices are payload and must survive sorting.
        // - After sorting by start we expect [10,15) idx 3, [15,18) idx 5, [20,30) idx 7.
        // - Span is therefore min start 10 and max end 30.
        let grouped = GroupedWindows::from_tuples(&[(20, 30, 7), (10, 15, 3), (15, 18, 5)], None)
            .expect("grouped test windows should be valid");

        let windows = grouped.windows_as_slice();

        assert_eq!(windows.len(), 3);
        assert_eq!(windows[0].into_tuple(), (10, 15, 3));
        assert_eq!(windows[1].into_tuple(), (15, 18, 5));
        assert_eq!(windows[2].into_tuple(), (20, 30, 7));
        assert_eq!(grouped.span_start(), 10);
        assert_eq!(grouped.span_end(), 30);
    }

    #[test]
    fn grouped_windows_sort_and_preserve_optional_strands() {
        // Arrange:
        // - Inputs are unsorted by start.
        // - Strand values are payload attached to rows, so they must move with their windows.
        // - After sorting by start, [10,15) keeps Forward, [15,18) keeps Unstranded,
        //   and [20,30) keeps Reverse.
        let grouped = GroupedWindows::new(
            vec![
                IndexedInterval::new(20, 30, 7).expect("grouped interval should be valid"),
                IndexedInterval::new(10, 15, 3).expect("grouped interval should be valid"),
                IndexedInterval::new(15, 18, 5).expect("grouped interval should be valid"),
            ],
            Some(vec![Strand::Reverse, Strand::Forward, Strand::Unstranded]),
        );

        let windows = grouped.windows_as_slice();
        let strands = grouped
            .strands
            .as_ref()
            .expect("strand metadata should be retained");

        assert_eq!(windows[0].into_tuple(), (10, 15, 3));
        assert_eq!(strands[0], Strand::Forward);
        assert_eq!(windows[1].into_tuple(), (15, 18, 5));
        assert_eq!(strands[1], Strand::Unstranded);
        assert_eq!(windows[2].into_tuple(), (20, 30, 7));
        assert_eq!(strands[2], Strand::Reverse);
    }

    #[test]
    fn grouped_windows_span_uses_max_end_not_last_sorted_end() {
        // Sorting by start yields [10,40), [20,25), [30,32). The last sorted window ends at 32,
        // but the collection span must use the true maximum end 40.
        let grouped = GroupedWindows::new(
            vec![
                IndexedInterval::new(20, 25, 0).expect("grouped interval should be valid"),
                IndexedInterval::new(10, 40, 1).expect("grouped interval should be valid"),
                IndexedInterval::new(30, 32, 2).expect("grouped interval should be valid"),
            ],
            None,
        );

        assert_eq!(grouped.span(), Span::new(10, 40).unwrap());
    }

    #[test]
    fn grouped_windows_empty_has_zero_span() {
        let grouped = GroupedWindows::from_sorted(Vec::new(), None);

        assert!(grouped.windows_as_slice().is_empty());
        assert_eq!(grouped.span_start(), 0);
        assert_eq!(grouped.span_end(), 0);
    }

    #[test]
    fn grouped_windows_from_tuples_rejects_invalid_interval() {
        let error = GroupedWindows::from_tuples(&[(10, 10, 3)], None)
            .expect_err("invalid grouped interval should fail");

        assert_eq!(
            error.to_string(),
            "interval end (10) must be greater than start (10)"
        );
    }

    #[test]
    fn scored_windows_sort_and_preserve_scores() {
        // Arrange:
        // - Sorting by start should reorder [20,30) and [10,15) into [10,15), [20,30).
        // - Score and original index are payload, so they must stay attached to their intervals.
        // - Span is the overall covered range [10,30).
        let scored = ScoredWindows::from_tuples(&[(20, 30, 7, 1.5), (10, 15, 3, 2.5)])
            .expect("scored test windows should be valid");

        let windows = scored.as_slice();

        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].into_tuple(), (10, 15, 3, 2.5));
        assert_eq!(windows[1].into_tuple(), (20, 30, 7, 1.5));
        assert_eq!(scored.span_start(), 10);
        assert_eq!(scored.span_end(), 30);
    }

    #[test]
    fn scored_windows_span_uses_max_end_not_last_sorted_end() {
        // Sorting by start yields [10,45), [20,25), [30,33). As with grouped windows, span_end
        // must be the global maximum end 45 rather than the last sorted end 33.
        let scored = ScoredWindows::new(vec![
            ScoredInterval::new(20, 25, 0, 0.5).expect("scored interval should be valid"),
            ScoredInterval::new(10, 45, 1, 1.5).expect("scored interval should be valid"),
            ScoredInterval::new(30, 33, 2, 2.5).expect("scored interval should be valid"),
        ]);

        assert_eq!(scored.span(), Span::new(10, 45).unwrap());
    }

    #[test]
    fn scored_windows_to_windows_drops_score_but_keeps_interval_and_index() {
        // Converting scored windows to plain windows should discard only the score field.
        // Interval bounds, original indices, and the collection span must remain unchanged.
        let scored = ScoredWindows::new(vec![
            ScoredInterval::new(5, 9, 11, 0.5).expect("scored interval should be valid"),
            ScoredInterval::new(10, 15, 12, 1.0).expect("scored interval should be valid"),
        ]);

        let plain = scored.to_windows();
        let windows = plain.as_slice();

        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].into_tuple(), (5, 9, 11));
        assert_eq!(windows[1].into_tuple(), (10, 15, 12));
        assert_eq!(plain.span(), Span::new(5, 15).unwrap());
    }

    #[test]
    fn windows_from_tuples_rejects_invalid_interval() {
        let error =
            Windows::from_tuples(&[(12, 12, 0)]).expect_err("invalid plain interval should fail");

        assert_eq!(
            error.to_string(),
            "interval end (12) must be greater than start (12)"
        );
    }

    #[test]
    fn scored_windows_from_tuples_rejects_invalid_interval() {
        let error = ScoredWindows::from_tuples(&[(20, 19, 7, 1.5)])
            .expect_err("invalid scored interval should fail");

        assert_eq!(
            error.to_string(),
            "interval end (19) must be greater than start (20)"
        );
    }
}
