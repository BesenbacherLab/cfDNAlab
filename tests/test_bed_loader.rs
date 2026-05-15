use anyhow::Result;
use cfdnalab::shared::bed::{
    GroupedBedStrandColumn, Strand, detect_header, line_looks_like_header,
    load_grouped_windows_from_bed, load_scored_windows_from_bed, load_windows_from_bed,
    write_group_idx_to_name_tsv,
};
use cfdnalab::shared::interval::{IndexedInterval, ScoredInterval};
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
fn should_error_when_column_5_looks_stranded_but_column_6_exists_without_strands() -> Result<()> {
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
fn should_treat_wide_grouped_bed_as_unstranded_when_no_strand_column_is_detected() -> Result<()> {
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

    let map = load_scored_windows_from_bed(bed.path(), Some(whitelist.as_slice()), None, Some(3))?;

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

#[test]
fn should_detect_header_after_comments_and_blank_lines() -> Result<()> {
    // The detector should skip comments and empty lines, so the first meaningful line here is the
    // literal header "chrom\tstart\tend".
    let bed = write_bed(&["# comment", "", "chrom\tstart\tend", "chr1\t0\t10"])?;

    assert!(detect_header(bed.path(), '\t')?);
    Ok(())
}

#[test]
fn should_detect_coordinate_lines_without_header() {
    // "chr1\t0\t10" has numeric coordinate columns and is therefore data, not a header.
    // The literal column names and comment lines should be treated as header-like.
    assert!(!line_looks_like_header("chr1\t0\t10", '\t'));
    assert!(line_looks_like_header("chrom\tstart\tend", '\t'));
    assert!(line_looks_like_header("# comment", '\t'));
}
