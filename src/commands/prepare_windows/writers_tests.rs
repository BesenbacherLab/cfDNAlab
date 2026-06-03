mod tests_writers {
    use crate::commands::prepare_windows::writers::{
        ChromTempWriter, ensure_temp_writer_for_chrom, finalize_temp_writers,
    };
    use fxhash::FxHashMap;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn ensure_temp_writer_creates_and_reuses_writer() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        {
            let writer = ensure_temp_writer_for_chrom("chr/1", dir.path(), &mut writers)?;
            writer.writer().write_all(b"chr1\t0\t5\n")?;
        }
        ensure_temp_writer_for_chrom("chr/1", dir.path(), &mut writers)?;
        assert_eq!(writers.len(), 1);
        let entries = finalize_temp_writers(&mut writers)?;
        let filename = entries[0]
            .1
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        assert_eq!(filename, "chrom.chrom-000000.bed.tmp");
        Ok(())
    }

    #[test]
    fn finalize_temp_writers_returns_empty_when_no_writers() -> anyhow::Result<()> {
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        let entries = finalize_temp_writers(&mut writers)?;
        assert!(entries.is_empty());
        Ok(())
    }

    #[test]
    fn finalize_temp_writers_flushes_and_clears() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        {
            let writer = ensure_temp_writer_for_chrom("chr1", dir.path(), &mut writers)?;
            writer.writer().write_all(b"chr1\t0\t5\n")?;
        }
        let entries = finalize_temp_writers(&mut writers)?;
        assert!(writers.is_empty());
        assert_eq!(entries.len(), 1);
        assert!(entries[0].1.exists());
        Ok(())
    }
}
