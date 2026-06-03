mod tests_filters {
    use super::super::filter_and_write_output;
    use crate::commands::prepare_windows::{
        config::PrepareConfig,
        labels::{AtomicLabelPart, LabelKey, LabelSchema},
    };
    use std::fs;
    use std::path::{Path, PathBuf};
    use tempfile::TempDir;

    fn write_temp_entries(
        base_dir: &Path,
        separator_name: &str,
        contents: &str,
    ) -> anyhow::Result<(PathBuf, Vec<(String, PathBuf)>)> {
        let stream_dir = base_dir.join(separator_name);
        fs::create_dir(&stream_dir)?;
        let temp_path = stream_dir.join("chr1.tmp");
        fs::write(&temp_path, contents)?;
        Ok((stream_dir, vec![("chr1".to_string(), temp_path)]))
    }

    fn input_label_schema() -> anyhow::Result<(LabelSchema, Vec<LabelKey>)> {
        let label_schema = LabelSchema::new(&[]).map_err(anyhow::Error::msg)?;
        let out_labels = vec![LabelKey::Atomic(AtomicLabelPart::Input)];
        Ok((label_schema, out_labels))
    }

    #[test]
    fn filter_and_write_output_replays_unfiltered_intermediate_rows() -> anyhow::Result<()> {
        let temp_dir = TempDir::new()?;
        let output = temp_dir.path().join("out.tsv");
        let (_stream_dir, entries) = write_temp_entries(
            temp_dir.path(),
            "stream",
            "chr1\t0\t5\tG||||\nchr1\t10\t15\tG||||\nchr1\t20\t25\tH||||\n",
        )?;
        let mut cfg = PrepareConfig::default();
        cfg.output = output.clone();
        cfg.sep = '\t';
        let (label_schema, out_labels) = input_label_schema()?;
        let chrom_order = vec!["chr1".to_string()];

        filter_and_write_output(
            &cfg,
            &entries,
            &label_schema,
            &out_labels,
            &[],
            &[],
            temp_dir.path(),
            &chrom_order,
        )?;

        assert_eq!(
            fs::read_to_string(output)?,
            "chr1\t0\t5\tG\nchr1\t10\t15\tG\nchr1\t20\t25\tH\n"
        );
        Ok(())
    }

    #[test]
    fn filter_and_write_output_uses_configured_separator_for_input_and_output()
    -> anyhow::Result<()> {
        let temp_dir = TempDir::new()?;
        let output = temp_dir.path().join("out.csv");
        let (_stream_dir, entries) =
            write_temp_entries(temp_dir.path(), "stream_csv", "chr1,0,5,G||||\n")?;
        let mut cfg = PrepareConfig::default();
        cfg.output = output.clone();
        cfg.sep = ',';
        let (label_schema, out_labels) = input_label_schema()?;
        let chrom_order = vec!["chr1".to_string()];

        filter_and_write_output(
            &cfg,
            &entries,
            &label_schema,
            &out_labels,
            &[],
            &[],
            temp_dir.path(),
            &chrom_order,
        )?;

        assert_eq!(fs::read_to_string(output)?, "chr1,0,5,G\n");
        Ok(())
    }
}
