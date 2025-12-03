use anyhow::Result;
use cfdnalab::shared::bed::load_windows_from_bed;
use flate2::{Compression, write::GzEncoder};
use std::io::Write;
use tempfile::NamedTempFile;

/// Write helper content to a temporary BED file used across tests.
fn write_bed(lines: &[&str]) -> Result<NamedTempFile> {
    let mut file = NamedTempFile::new()?;
    for line in lines {
        writeln!(file, "{}", line)?;
    }
    Ok(file)
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
    assert_eq!(chr1.as_slice(), &[(0, 10, 0), (20, 30, 2)]);

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
    assert_eq!(chr1.as_slice(), &[(10, 25, 1)]);
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
    assert_eq!(chr1.as_slice(), &[(0, 5, 0), (10, 15, 1)]);
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
    assert_eq!(chr2.as_slice(), &[(4, 8, 1), (8, 12, 2)]);

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
