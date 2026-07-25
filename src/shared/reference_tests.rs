use super::*;
use crate::testing::reference::twobit_from_sequences;

#[test]
fn staged_references_are_distinct_copies_with_identical_bytes() -> anyhow::Result<()> {
    let reference = twobit_from_sequences(
        "staged_reference",
        vec![("chr1".to_string(), "ACGTACGTACGT".to_string())],
    )?;
    let first_work_dir = tempfile::TempDir::new()?;
    let second_work_dir = tempfile::TempDir::new()?;

    let first_copy = stage_reference_2bit(reference.path(), first_work_dir.path())?;
    let second_copy = stage_reference_2bit(reference.path(), second_work_dir.path())?;

    assert_eq!(fs::read(reference.path())?, fs::read(&first_copy)?);
    assert_eq!(fs::read(&first_copy)?, fs::read(&second_copy)?);
    assert_ne!(first_copy, second_copy);
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        assert_ne!(
            fs::metadata(&first_copy)?.ino(),
            fs::metadata(&second_copy)?.ino()
        );
    }
    Ok(())
}

#[test]
fn reusable_reader_supports_repeated_full_and_ranged_reads() -> anyhow::Result<()> {
    let reference = twobit_from_sequences(
        "reusable_reference_reader",
        vec![
            ("chr1".to_string(), "ACGTACGT".to_string()),
            ("chr2".to_string(), "TTGCAACC".to_string()),
        ],
    )?;
    let mut reader = ReferenceReader::open(reference.path())?;

    assert_eq!(reader.read_seq("chr1")?, b"ACGTACGT");
    assert_eq!(reader.read_seq_in_range("chr1", 2..6)?, b"GTAC");
    assert_eq!(reader.read_seq("chr2")?, b"TTGCAACC");
    assert_eq!(reader.read_seq_in_range("chr1", 0..4)?, b"ACGT");
    Ok(())
}
