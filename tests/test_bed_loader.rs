use anyhow::Result;
use cfdnalab::shared::bed::load_windows_from_bed;
use tempfile::NamedTempFile;

/// Write helper content to a temporary BED file used across tests.
fn write_bed(lines: &[&str]) -> Result<NamedTempFile> {
    let mut file = NamedTempFile::new()?;
    use std::io::Write;
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
    let map = load_windows_from_bed(bed.path(), Some(whitelist.as_slice()), None)?;

    // Assert
    let chr1 = map.get("chr1").expect("chr1 missing");
    assert_eq!(chr1.as_slice(), &[(0, 10, 0), (20, 30, 1)]);

    let empty = load_windows_from_bed(bed.path(), Some(["chr3".to_string()].as_slice()), None)?;
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
    let map = load_windows_from_bed(bed.path(), None, Some(&keep_large))?;

    // Assert
    let chr1 = map.get("chr1").expect("chr1 missing");
    assert_eq!(chr1.as_slice(), &[(10, 25, 0)]);
    Ok(())
}
