mod tests_chunk {

    use crate::commands::prepare_windows::{
        chunk::{flush_chromosome, process_and_write_chunk},
        config::{DedupKeep, DistancePolicy, MergeLabel, MergeScope, PrepareConfig},
        intermediate::parse_intermediate_line,
        labels::{AtomicLabelPart, LabelKey, LabelSchema, LabelTuple},
        prepare_windows::Window,
        writers::{ChromTempWriter, finalize_temp_writers},
    };
    use fxhash::FxHashMap;
    use std::fs;
    use std::sync::Arc;
    use tempfile::TempDir;

    fn win(chrom: &str, start: u32, end: u32, group: &str) -> Window {
        Window::from_bounds(
            Arc::<str>::from(chrom.to_string()),
            start,
            end,
            start,
            end,
            vec![LabelTuple::new(group.to_string())],
            group.to_string(),
            None,
        )
        .expect("test window should be valid")
    }

    fn label_schema() -> LabelSchema {
        LabelSchema::new(&[]).expect("label schema")
    }

    fn merge_key() -> LabelKey {
        LabelKey::Atomic(AtomicLabelPart::Input)
    }

    fn out_labels() -> Vec<LabelKey> {
        vec![LabelKey::Atomic(AtomicLabelPart::Input)]
    }

    fn make_config() -> PrepareConfig {
        let mut cfg = PrepareConfig::default();
        cfg.deduplicate = DedupKeep::None;
        cfg.min_distance_within_group = None;
        cfg.distance_policy = DistancePolicy::KeepFirst;
        cfg.merge_scope = MergeScope::None;
        cfg.merge_gap = None;
        cfg.merge_label = MergeLabel::Join;
        cfg.sep = '\t';
        cfg
    }

    #[test]
    fn process_and_write_chunk_writes_safe_prefix() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let cfg = make_config();
        let mut carryover = Vec::new();
        let mut batch = vec![win("chr1", 0, 5, "g"), win("chr1", 10, 15, "g")];
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        let schema = label_schema();
        let key = merge_key();
        let labels = out_labels();
        let mut near_index = None;
        process_and_write_chunk(
            "chr1",
            &mut carryover,
            &mut batch,
            &mut writers,
            dir.path(),
            None,
            0,
            None,
            &cfg,
            &mut near_index,
            None,
            &schema,
            &key,
            &labels,
        )?;
        assert!(carryover.is_empty());
        let entries = finalize_temp_writers(&mut writers)?;
        let contents = fs::read_to_string(&entries[0].1)?;
        assert!(contents.contains("chr1\t0\t5"));
        assert!(contents.contains("chr1\t10\t15"));
        Ok(())
    }

    #[test]
    fn process_and_write_chunk_retains_tail_when_margin_present() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let mut cfg = make_config();
        cfg.min_distance_within_group = Some(5);
        cfg.merge_scope = MergeScope::Within;
        let mut carryover = Vec::new();
        let mut batch = vec![win("chr1", 0, 5, "g"), win("chr1", 3, 8, "g")];
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        let schema = label_schema();
        let key = merge_key();
        let labels = out_labels();
        let mut near_index = None;
        process_and_write_chunk(
            "chr1",
            &mut carryover,
            &mut batch,
            &mut writers,
            dir.path(),
            None,
            0,
            None,
            &cfg,
            &mut near_index,
            None,
            &schema,
            &key,
            &labels,
        )?;
        assert_eq!(carryover.len(), 1);
        let entries = finalize_temp_writers(&mut writers)?;
        let contents = fs::read_to_string(&entries[0].1)?;
        assert!(contents.is_empty());
        Ok(())
    }

    #[test]
    fn flush_chromosome_writes_remaining_tail() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let mut cfg = make_config();
        cfg.min_distance_within_group = Some(5);
        let mut carryover = vec![win("chr1", 0, 5, "g")];
        let mut batch = vec![win("chr1", 5, 9, "g")];
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        let schema = label_schema();
        let key = merge_key();
        let labels = out_labels();
        let mut near_index = None;
        flush_chromosome(
            "chr1",
            &mut carryover,
            &mut batch,
            &mut writers,
            dir.path(),
            None,
            0,
            None,
            &cfg,
            &mut near_index,
            None,
            &schema,
            &key,
            &labels,
        )?;
        assert!(carryover.is_empty());
        let entries = finalize_temp_writers(&mut writers)?;
        let contents = fs::read_to_string(&entries[0].1)?;
        assert!(contents.contains("chr1\t0\t5"));
        assert!(contents.contains("chr1\t5\t9"));
        Ok(())
    }

    #[test]
    fn process_and_write_chunk_applies_deduplication() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let mut cfg = make_config();
        cfg.deduplicate = DedupKeep::KeepFirst;
        let mut carryover = Vec::new();
        let mut batch = vec![win("chr1", 0, 5, "g"), win("chr1", 0, 5, "g")];
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        let schema = label_schema();
        let key = merge_key();
        let labels = out_labels();
        let mut near_index = None;
        process_and_write_chunk(
            "chr1",
            &mut carryover,
            &mut batch,
            &mut writers,
            dir.path(),
            None,
            0,
            None,
            &cfg,
            &mut near_index,
            None,
            &schema,
            &key,
            &labels,
        )?;
        let entries = finalize_temp_writers(&mut writers)?;
        let contents = fs::read_to_string(&entries[0].1)?;
        assert_eq!(contents.trim().lines().count(), 1);
        Ok(())
    }

    #[test]
    fn flush_chromosome_is_noop_when_empty() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let cfg = make_config();
        let mut carryover = Vec::new();
        let mut batch = Vec::new();
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        let schema = label_schema();
        let key = merge_key();
        let labels = out_labels();
        let mut near_index = None;
        flush_chromosome(
            "chr1",
            &mut carryover,
            &mut batch,
            &mut writers,
            dir.path(),
            None,
            0,
            None,
            &cfg,
            &mut near_index,
            None,
            &schema,
            &key,
            &labels,
        )?;
        assert!(writers.is_empty());
        Ok(())
    }

    #[test]
    fn chunking_merges_across_scope_over_chunk_boundary() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let mut cfg = make_config();
        cfg.merge_scope = MergeScope::Across;
        cfg.merge_gap = Some(2);
        let mut carryover = Vec::new();
        let mut batch = vec![win("chr1", 0, 5, "g1"), win("chr1", 7, 10, "g2")];
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        let schema = label_schema();
        let key = merge_key();
        let labels = out_labels();
        let mut near_index = None;

        process_and_write_chunk(
            "chr1",
            &mut carryover,
            &mut batch,
            &mut writers,
            dir.path(),
            None,
            0,
            None,
            &cfg,
            &mut near_index,
            None,
            &schema,
            &key,
            &labels,
        )?;

        assert_eq!(carryover.len(), 1); // retained tail for next chunk

        // Flush remaining tail and ensure merged output is written
        flush_chromosome(
            "chr1",
            &mut carryover,
            &mut Vec::new(),
            &mut writers,
            dir.path(),
            None,
            0,
            None,
            &cfg,
            &mut near_index,
            None,
            &schema,
            &key,
            &labels,
        )?;

        let entries = finalize_temp_writers(&mut writers)?;
        let contents = fs::read_to_string(&entries[0].1)?;
        let line = contents.lines().next().expect("intermediate line");
        let parsed = parse_intermediate_line(line, cfg.sep)?;
        let inputs: Vec<&str> = parsed
            .label_tuples
            .iter()
            .map(|tuple| tuple.input.as_str())
            .collect();
        // Tuples are stored separately in intermediate files
        assert_eq!(parsed.chrom, "chr1");
        assert_eq!(parsed.start(), 0);
        assert_eq!(parsed.end(), 10);
        assert_eq!(inputs, vec!["g1", "g2"]);
        Ok(())
    }
}
