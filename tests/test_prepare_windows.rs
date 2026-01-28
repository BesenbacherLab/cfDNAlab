#![cfg(feature = "cmd_prepare_windows")]

mod tests_prepare_windows_pipeline {

    use anyhow::Result;
    use cfdnalab::commands::prepare_windows::config::{
        CoordinateSet, DedupKeep, DistSign, DistancePolicy, HeaderMode, MergeLabel, MergeScope,
        NearDirection, NearEdge, NearTiePolicy, OobPolicy, PrepareConfig,
    };
    use cfdnalab::commands::prepare_windows::prepare_windows::run;
    use flate2::{Compression, write::GzEncoder};
    use std::fs;
    use std::io::{BufWriter, Read, Write};
    use tempfile::{NamedTempFile, TempDir};
    use zstd::Decoder as ZstdDecoder;

    fn write_temp_file(dir: &TempDir, name: &str, lines: &[&str]) -> Result<std::path::PathBuf> {
        let path = dir.path().join(name);
        fs::write(&path, lines.join("\n") + "\n")?;
        Ok(path)
    }

    fn write_chrom_sizes(dir: &TempDir, sizes: &[(&str, u32)]) -> Result<std::path::PathBuf> {
        let path = dir.path().join("chrom.sizes");
        let content = sizes
            .iter()
            .map(|(chrom, size)| format!("{chrom}\t{size}"))
            .collect::<Vec<_>>()
            .join("\n");
        fs::write(&path, content + "\n")?;
        Ok(path)
    }

    fn run_pipeline(cfg: &PrepareConfig) -> Result<Vec<String>> {
        run(cfg)?;
        let output_path = cfg.output.clone();
        let contents = fs::read_to_string(output_path)?;
        Ok(contents
            .lines()
            .filter(|line| !line.is_empty())
            .map(|line| line.to_string())
            .collect())
    }

    #[test]
    fn should_write_windows_verbatim_when_no_transformations() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "windows.tsv", &["chr1\t0\t10", "chr1\t10\t20"])?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output.clone();
        cfg.header = HeaderMode::Absent;
        cfg.oob = OobPolicy::Allow;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(lines, vec!["chr1\t0\t10\t", "chr1\t10\t20\t"]);
        Ok(())
    }

    #[test]
    fn should_drop_blacklisted_windows_and_label_with_near_distance() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &[
                "chr1\t0\t10\tG1",
                "chr1\t15\t25\tG1",
                "chr1\t40\t50\tG2",
                "chr2\t5\t15\tG2",
            ],
        )?;
        let blacklist = write_temp_file(&tmpdir, "blacklist.bed", &["chr1\t8\t18"])?;
        let near = write_temp_file(
            &tmpdir,
            "near.bed",
            &["chr1\t35\t45\tSiteA", "chr2\t0\t10\tSiteB"],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.blacklist = Some(vec![blacklist]);
        cfg.blacklist_strategy = cfdnalab::shared::blacklist::BlacklistStrategy::Any;
        cfg.blacklist_halo = 0;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.near_edge = NearEdge::Nearest;
        cfg.near_direction = NearDirection::Both;
        cfg.distance_sign = DistSign::Absolute;
        cfg.distance_bins = Some(vec!["prox:<10".to_string(), "far:>=10".to_string()]);
        cfg.group_cols = vec!["3".to_string()];
        cfg.out_labels = vec![
            "input".to_string(),
            "near-side".to_string(),
            "near-name".to_string(),
            "bin".to_string(),
        ];
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(
            lines,
            vec![
                "chr1\t40\t50\tG2\t=\tSiteA\tprox".to_string(),
                "chr2\t5\t15\tG2\t=\tSiteB\tprox".to_string(),
            ],
        );
        Ok(())
    }

    #[test]
    fn should_merge_then_apply_blacklist() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t10\t14\tG", "chr1\t18\t22\tG", "chr1\t40\t44\tH"],
        )?;
        let blacklist = write_temp_file(&tmpdir, "blacklist.bed", &["chr1\t15\t17"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.blacklist = Some(vec![blacklist]);
        cfg.blacklist_strategy = cfdnalab::shared::blacklist::BlacklistStrategy::Any;
        cfg.group_cols = vec!["3".to_string()];
        cfg.merge_scope = MergeScope::Within;
        cfg.merge_gap = Some(4);
        cfg.merge_label = MergeLabel::Join;
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;
        assert_eq!(lines, vec!["chr1\t40\t44\tH".to_string()]);
        Ok(())
    }

    #[test]
    fn should_resize_before_merging() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t10\t14\tA", "chr1\t20\t24\tA", "chr1\t40\t44\tA"],
        )?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.flank = Some(vec![3, 3]);
        cfg.oob = OobPolicy::Allow;
        cfg.merge_scope = MergeScope::Within;
        cfg.merge_gap = Some(0);
        cfg.merge_label = MergeLabel::Join;
        cfg.merge_on = CoordinateSet::Resized;
        cfg.chrom_sizes = Some(chrom_sizes);

        let lines = run_pipeline(&cfg)?;
        assert_eq!(
            lines,
            vec!["chr1\t7\t27\tA".to_string(), "chr1\t37\t47\tA".to_string(),],
        );
        Ok(())
    }

    #[test]
    fn should_resize_after_merging_when_merge_on_original() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t10\t14\tA", "chr1\t20\t24\tA", "chr1\t40\t44\tA"],
        )?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.resize = Some(5);
        cfg.oob = OobPolicy::Allow;
        cfg.merge_scope = MergeScope::Within;
        cfg.merge_gap = Some(6);
        cfg.merge_label = MergeLabel::Join;
        cfg.merge_on = CoordinateSet::Original;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        // Merge first two windows because original gap is 6 and merge_gap is 6 (10->24)
        // Midpoint is start + (end-start)/2, so 10-24 -> midpoint 17, resize 5 -> 15-20
        assert_eq!(
            lines,
            vec!["chr1\t15\t20\tA".to_string(), "chr1\t40\t45\tA".to_string(),],
        );
        Ok(())
    }

    #[test]
    fn should_assign_three_distance_bins() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t6\t8", "chr1\t52\t54", "chr1\t150\t152"],
        )?;
        let near = write_temp_file(
            &tmpdir,
            "near.tsv",
            &["chr1\t10\t12\tA", "chr1\t40\t42\tB", "chr1\t80\t82\tC"],
        )?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_edge = NearEdge::Nearest;
        cfg.distance_sign = DistSign::Absolute;
        cfg.distance_bins = Some(vec![
            "prox:<5".to_string(),
            "mid:5-15".to_string(),
            "far:>15".to_string(),
        ]);
        cfg.out_labels = vec!["near-side".to_string(), "bin".to_string()];
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;

        assert_eq!(
            lines,
            vec![
                "chr1\t6\t8\t+\tprox".to_string(),
                "chr1\t52\t54\t-\tmid".to_string(),
                "chr1\t150\t152\t-\tfar".to_string(),
            ],
        );
        Ok(())
    }

    #[test]
    fn should_use_resized_coordinates_for_distance_binning() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t20\t30"])?;
        let near = write_temp_file(&tmpdir, "near.tsv", &["chr1\t0\t10\tSITE"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.near_edge = NearEdge::Nearest;
        cfg.near_direction = NearDirection::Both;
        cfg.distance_sign = DistSign::Absolute;
        cfg.distance_bins = Some(vec!["near:<7".to_string(), "far:>=7".to_string()]);
        cfg.distance_from = CoordinateSet::Resized;
        cfg.flank = Some(vec![5, 5]);
        cfg.out_labels = vec![
            "near-side".to_string(),
            "near-name".to_string(),
            "bin".to_string(),
        ];
        cfg.oob = OobPolicy::Allow;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        // Resized 15..35 yields a 5 bp distance bin and output keeps resized coordinates
        assert_eq!(lines, vec!["chr1\t15\t35\t-\tSITE\tnear".to_string()]);
        Ok(())
    }

    #[test]
    fn should_use_original_coordinates_for_distance_binning() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t20\t30"])?;
        let near = write_temp_file(&tmpdir, "near.tsv", &["chr1\t0\t10\tSITE"])?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.near_edge = NearEdge::Nearest;
        cfg.near_direction = NearDirection::Both;
        cfg.distance_sign = DistSign::Absolute;
        cfg.distance_bins = Some(vec!["near:<7".to_string(), "far:>=7".to_string()]);
        cfg.distance_from = CoordinateSet::Original;
        cfg.flank = Some(vec![5, 5]);
        cfg.out_labels = vec![
            "near-side".to_string(),
            "near-name".to_string(),
            "bin".to_string(),
        ];
        cfg.oob = OobPolicy::Allow;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        // Original 20..30 yields a 10 bp distance bin while output keeps resized coordinates
        assert_eq!(lines, vec!["chr1\t15\t35\t-\tSITE\tfar".to_string()]);
        Ok(())
    }

    #[test]
    fn should_annotate_ties_with_directional_groups() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t10\t20"])?;
        let near = write_temp_file(&tmpdir, "near.tsv", &["chr1\t0\t5\tUP", "chr1\t25\t30\tDN"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.near_edge = NearEdge::Nearest;
        cfg.near_direction = NearDirection::Both;
        cfg.distance_sign = DistSign::Signed;
        cfg.out_labels = vec!["near-side".to_string(), "near-name".to_string()];
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;
        assert_eq!(lines, vec!["chr1\t10\t20\t+,-\tDN,UP".to_string()]);
        Ok(())
    }

    #[test]
    fn should_drop_ties_when_configured() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t10\t20"])?;
        let near = write_temp_file(&tmpdir, "near.tsv", &["chr1\t0\t5\tUP", "chr1\t25\t30\tDN"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_edge = NearEdge::Nearest;
        cfg.near_direction = NearDirection::Both;
        cfg.distance_sign = DistSign::Signed;
        cfg.near_ties = NearTiePolicy::Drop;
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;
        assert!(lines.is_empty());
        Ok(())
    }

    #[test]
    fn should_add_direction_prefix_for_unique_hits() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t0\t10", "chr1\t150\t160"])?;
        let near = write_temp_file(&tmpdir, "near.tsv", &["chr1\t100\t110\tTARG"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.near_edge = NearEdge::Nearest;
        cfg.distance_sign = DistSign::Absolute; // Direction prefix should still appear
        cfg.out_labels = vec!["near-side".to_string(), "near-name".to_string()];
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;

        assert_eq!(
            lines,
            vec!["chr1\t0\t10\t+\tTARG", "chr1\t150\t160\t-\tTARG"],
        );
        Ok(())
    }

    #[test]
    fn should_add_direction_prefix_without_near_groups() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t0\t10", "chr1\t150\t160"])?;
        let near = write_temp_file(&tmpdir, "near.tsv", &["chr1\t100\t110\tTARG"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_edge = NearEdge::Nearest;
        cfg.distance_sign = DistSign::Absolute; // Direction prefix should still appear
        cfg.out_labels = vec!["near-side".to_string()];
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(
            lines,
            vec![
                "chr1\t0\t10\t+".to_string(),
                "chr1\t150\t160\t-".to_string()
            ]
        );
        Ok(())
    }

    #[test]
    fn should_exclude_windows_by_atomic_label() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t0\t10\tA", "chr1\t10\t20\tB", "chr1\t20\t30\tA"],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.exclude_labels = vec!["input=B".to_string()];
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(
            lines,
            vec!["chr1\t0\t10\tA".to_string(), "chr1\t20\t30\tA".to_string()],
        );
        Ok(())
    }

    #[test]
    fn should_exclude_windows_by_composition_label() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t0\t10\tA", "chr1\t40\t50\tB"])?;
        let near = write_temp_file(
            &tmpdir,
            "near.tsv",
            &["chr1\t0\t5\t+\tSite1", "chr1\t40\t45\t+\tSite2"],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.compose = vec!["core=input,bin".parse().map_err(anyhow::Error::msg)?];
        cfg.out_labels = vec!["input".to_string(), "core".to_string()];
        cfg.exclude_labels = vec!["core=A.prox".to_string()];
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_edge = NearEdge::Nearest;
        cfg.near_direction = NearDirection::Both;
        cfg.distance_sign = DistSign::Absolute;
        cfg.distance_bins = Some(vec!["prox:<5".to_string()]);
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(lines, vec!["chr1\t40\t50\tB\tB.prox".to_string()]);
        Ok(())
    }

    #[test]
    fn should_filter_by_min_per_for_input_groups() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t0\t10\tA", "chr1\t10\t20\tA", "chr1\t20\t30\tB"],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.out_labels = vec!["input".to_string()];
        cfg.min_per = vec!["input=2".to_string()];
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(
            lines,
            vec!["chr1\t0\t10\tA".to_string(), "chr1\t10\t20\tA".to_string()],
        );
        Ok(())
    }

    #[test]
    fn should_tag_clusters_across_groups() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &[
                "chr1\t0\t10\tA",
                "chr1\t0\t10\tB",
                "chr1\t0\t10\tC",
                "chr1\t20\t30\tD",
            ],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.out_labels = vec!["input".to_string(), "cluster".to_string()];
        cfg.cluster_min_overlaps = Some(2);
        cfg.cluster_before_min_distance = false;
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(
            lines,
            vec![
                "chr1\t0\t10\tA\tcluster".to_string(),
                "chr1\t0\t10\tB\tcluster".to_string(),
                "chr1\t0\t10\tC\tcluster".to_string(),
                "chr1\t20\t30\tD\tnone".to_string(),
            ],
        );
        Ok(())
    }

    #[test]
    fn should_cluster_before_min_distance_when_configured() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t0\t10\tA", "chr1\t0\t10\tA"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.out_labels = vec!["input".to_string(), "cluster".to_string()];
        cfg.cluster_min_overlaps = Some(2);
        cfg.cluster_before_min_distance = true;
        cfg.min_distance_within_group = Some(1);
        cfg.distance_policy = DistancePolicy::KeepFirst;
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(lines, vec!["chr1\t0\t10\tA\tcluster".to_string()]);
        Ok(())
    }

    #[test]
    fn should_cluster_after_min_distance_by_default() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t0\t10\tA", "chr1\t0\t10\tA"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.group_cols = vec!["3".to_string()];
        cfg.out_labels = vec!["input".to_string(), "cluster".to_string()];
        cfg.cluster_min_overlaps = Some(2);
        cfg.cluster_before_min_distance = false;
        cfg.min_distance_within_group = Some(1);
        cfg.distance_policy = DistancePolicy::KeepFirst;
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(lines, vec!["chr1\t0\t10\tA\tnone".to_string()]);
        Ok(())
    }

    #[test]
    fn should_resize_and_deduplicate_windows_with_spacing() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t10\t20", "chr1\t10\t20", "chr1\t20\t22"],
        )?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.resize = Some(8);
        cfg.deduplicate = DedupKeep::KeepFirst;
        cfg.min_distance_within_group = Some(4);
        cfg.distance_policy = DistancePolicy::KeepFirst;
        cfg.oob = OobPolicy::Allow;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        // Resizing centers on the true midpoint (even->even has a unique placement).
        assert_eq!(lines, vec!["chr1\t11\t19\t", "chr1\t17\t25\t"]);
        Ok(())
    }

    #[test]
    fn should_flank_windows_and_clip_to_allowed_bounds() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t5\t10"])?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.flank = Some(vec![3, 7]);
        cfg.oob = OobPolicy::Allow;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(lines, vec!["chr1\t2\t17\t"]);
        Ok(())
    }

    #[test]
    fn should_merge_within_group() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &[
                "chr1\t0\t5\tG",
                "chr1\t4\t8\tG",
                "chr1\t50\t55\tH",
                "chr1\t60\t65\tH",
            ],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.merge_scope = MergeScope::Within;
        cfg.merge_gap = Some(2);
        cfg.merge_label = MergeLabel::Join;
        cfg.oob = OobPolicy::Allow;
        cfg.group_cols = vec!["3".to_string()];

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(
            lines,
            vec!["chr1\t0\t8\tG", "chr1\t50\t55\tH", "chr1\t60\t65\tH"]
        );
        Ok(())
    }

    #[test]
    fn should_execute_full_pipeline_end_to_end() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &[
                "chr1\t10\t15\tA",
                "chr1\t16\t21\tA",
                "chr1\t40\t45\tB",
                "chr2\t0\t5\tC",
            ],
        )?;
        let blacklist = write_temp_file(&tmpdir, "mask.bed", &["chr1\t38\t50"])?;
        let near = write_temp_file(&tmpdir, "near.tsv", &["chr1\t8\t12\tN1", "chr2\t0\t10\tN2"])?;
        let output = tmpdir.path().join("out.tsv");
        let chrom_sizes = write_chrom_sizes(&tmpdir, &[("chr1", 1000), ("chr2", 1000)])?;

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.blacklist = Some(vec![blacklist]);
        cfg.blacklist_strategy = cfdnalab::shared::blacklist::BlacklistStrategy::Any;
        cfg.blacklist_halo = 2;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.near_edge = NearEdge::Nearest;
        cfg.distance_sign = DistSign::Signed;
        cfg.distance_bins = Some(vec!["close:<5".to_string(), "far:>=5".to_string()]);
        cfg.resize = Some(10);
        cfg.min_distance_within_group = Some(5);
        cfg.distance_policy = DistancePolicy::KeepLongest;
        cfg.merge_scope = MergeScope::Across;
        cfg.merge_gap = Some(3);
        cfg.merge_label = MergeLabel::Join;
        cfg.oob = OobPolicy::Allow;
        cfg.group_cols = vec!["3".to_string()];
        cfg.out_labels = vec![
            "input".to_string(),
            "near-side".to_string(),
            "near-name".to_string(),
            "bin".to_string(),
        ];
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        let expected_left = vec![
            "chr1\t7\t24\tA\t=\tN1\tclose".to_string(),
            "chr2\t0\t8\tC\t=\tN2\tclose".to_string(),
        ];
        // When resize parity is ambiguous, placement can shift by one bp depending on the
        // deterministic hash choice (left vs right). Accept either centred placement.
        let expected_right = vec![
            "chr1\t8\t24\tA\t=\tN1\tclose".to_string(),
            "chr2\t0\t8\tC\t=\tN2\tclose".to_string(),
        ];
        assert!(
            lines == expected_left || lines == expected_right,
            "unexpected output: {:?}",
            lines
        );
        Ok(())
    }

    #[test]
    fn should_error_when_chromosome_reappears_out_of_order() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(
            &tmpdir,
            "input.tsv",
            &["chr1\t0\t5", "chr2\t0\t5", "chr1\t10\t15"],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.oob = OobPolicy::Allow;

        // Act
        let err = run(&cfg).unwrap_err();

        // Assert
        assert!(
            format!("{err}").contains("chromosome 'chr1' appears after it was already processed")
        );
        Ok(())
    }

    #[test]
    fn should_error_when_start_coordinate_decreases() -> Result<()> {
        // Arrange
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t5\t10", "chr1\t3\t8"])?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.oob = OobPolicy::Allow;

        // Act
        let err = run(&cfg).unwrap_err();

        // Assert
        assert!(format!("{err}").contains("has start 3 before previous 5"));
        Ok(())
    }

    #[test]
    fn should_accept_gz_input_and_emit_zst_output() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let gz_path = tmpdir.path().join("input.tsv.gz");
        {
            let file = std::fs::File::create(&gz_path)?;
            let buf = BufWriter::new(file);
            let mut encoder = GzEncoder::new(buf, Compression::default());
            writeln!(encoder, "chr1\t0\t10")?;
            writeln!(encoder, "chr1\t10\t20")?;
            let mut buf = encoder.finish()?;
            buf.flush()?;
        }

        let output = tmpdir.path().join("out.tsv.zst");

        let mut cfg = PrepareConfig::default();
        cfg.input = gz_path;
        cfg.output = output.clone();
        cfg.header = HeaderMode::Absent;
        cfg.oob = OobPolicy::Allow;

        run(&cfg)?;

        let file = std::fs::File::open(&output)?;
        let mut decoder = ZstdDecoder::new(file)?;
        let mut text = String::new();
        decoder.read_to_string(&mut text)?;
        let lines: Vec<_> = text.lines().filter(|line| !line.is_empty()).collect();

        assert_eq!(lines, vec!["chr1\t0\t10\t", "chr1\t10\t20\t"]);
        Ok(())
    }
}

mod tests_postprocess {
    use cfdnalab::commands::prepare_windows::{
        config::{CoordinateSet, DedupKeep, DistancePolicy, MergeScope},
        labels::LabelTuple,
        postprocess::{
            deduplicate_identical, enforce_min_distance_within_group, partition_safe_and_tail,
        },
        prepare_windows::Window,
    };
    use std::sync::Arc;

    fn win(chrom: &str, start: u32, end: u32, group: &str, score: Option<f32>) -> Window {
        Window {
            chrom: Arc::<str>::from(chrom.to_string()),
            original_start: start,
            original_end: end,
            resized_start: start,
            resized_end: end,
            label_tuples: vec![LabelTuple::new(group.to_string())],
            group_key: group.to_string(),
            score,
            merged: false,
        }
    }

    fn snapshot(windows: &[Window]) -> Vec<(String, u32, u32, String, Option<f32>)> {
        windows
            .iter()
            .map(|w| {
                (
                    w.chrom.as_ref().to_string(),
                    w.resized_start,
                    w.resized_end,
                    w.group_key.clone(),
                    w.score,
                )
            })
            .collect()
    }

    #[test]
    fn dedup_none_keeps_all_windows() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(2.0)),
        ];
        let result = deduplicate_identical(
            windows.clone(),
            DedupKeep::None,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(snapshot(&result), snapshot(&windows));
    }

    #[test]
    fn dedup_keep_first_prefers_first_duplicate() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(5.0)),
        ];
        let result =
            deduplicate_identical(windows, DedupKeep::KeepFirst, true, CoordinateSet::Resized);
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", Some(1.0))])
        );
    }

    #[test]
    fn dedup_keep_highest_score_uses_scores_when_available() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(5.0)),
            win("chr1", 10, 20, "g1", Some(2.5)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", Some(5.0))])
        );
    }

    #[test]
    fn dedup_keep_highest_score_falls_back_without_scores() {
        let windows = vec![
            win("chr1", 10, 20, "g1", None),
            win("chr1", 10, 20, "g1", None),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            false,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", None)])
        );
    }

    #[test]
    fn dedup_keep_lowest_score_picks_smallest_score() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(3.0)),
            win("chr1", 10, 20, "g1", Some(1.5)),
            win("chr1", 10, 20, "g1", Some(4.0)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepLowestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", Some(1.5))])
        );
    }

    #[test]
    fn dedup_keep_lowest_score_falls_back_without_scores() {
        let windows = vec![
            win("chr1", 10, 20, "g1", None),
            win("chr1", 10, 20, "g1", None),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepLowestScore,
            false,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", None)])
        );
    }

    #[test]
    fn dedup_does_not_touch_unique_windows() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 30, 40, "g1", Some(2.0)),
            win("chr2", 5, 15, "", None),
        ];
        let result = deduplicate_identical(
            windows.clone(),
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(snapshot(&result), snapshot(&windows));
    }

    #[test]
    fn dedup_handles_multiple_duplicate_groups() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(3.0)),
            win("chr1", 30, 40, "g2", Some(5.0)),
            win("chr1", 30, 40, "g2", Some(2.0)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[
                win("chr1", 10, 20, "g1", Some(3.0)),
                win("chr1", 30, 40, "g2", Some(5.0)),
            ])
        );
    }

    #[test]
    fn dedup_keep_highest_score_prefers_non_none_scores() {
        let windows = vec![
            win("chr1", 0, 5, "g", None),
            win("chr1", 0, 5, "g", Some(1.0)),
            win("chr1", 0, 5, "g", Some(2.0)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 0, 5, "g", Some(2.0))])
        );
    }

    #[test]
    fn enforce_min_distance_within_group_keep_first() {
        let windows = vec![
            win("chr1", 0, 10, "g", Some(1.0)),
            win("chr1", 4, 12, "g", Some(2.0)),
            win("chr1", 20, 30, "g", Some(3.0)),
        ];
        let result = enforce_min_distance_within_group(
            windows,
            Some(5),
            DistancePolicy::KeepFirst,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[
                win("chr1", 0, 10, "g", Some(1.0)),
                win("chr1", 20, 30, "g", Some(3.0))
            ])
        );
    }

    #[test]
    fn enforce_min_distance_within_group_keep_highest_score() {
        let windows = vec![
            win("chr1", 0, 10, "g", Some(1.0)),
            win("chr1", 4, 12, "g", Some(5.0)),
            win("chr1", 40, 50, "g", Some(2.0)),
        ];
        let result = enforce_min_distance_within_group(
            windows,
            Some(8),
            DistancePolicy::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[
                win("chr1", 4, 12, "g", Some(5.0)),
                win("chr1", 40, 50, "g", Some(2.0))
            ])
        );
    }

    #[test]
    fn enforce_min_distance_within_group_keep_lowest_score_without_scores() {
        let windows = vec![
            win("chr1", 0, 5, "g", None),
            win("chr1", 3, 9, "g", None),
            win("chr1", 20, 25, "g", None),
        ];
        let result = enforce_min_distance_within_group(
            windows,
            Some(4),
            DistancePolicy::KeepLowestScore,
            false,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 0, 5, "g", None), win("chr1", 20, 25, "g", None)])
        );
    }

    #[test]
    fn partition_safe_and_tail_without_margin_writes_all() {
        let windows = vec![
            win("chr1", 0, 10, "g", None),
            win("chr1", 20, 30, "g", None),
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::Within,
            None,
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        assert_eq!(
            snapshot(&safe),
            snapshot(&[
                win("chr1", 0, 10, "g", None),
                win("chr1", 20, 30, "g", None)
            ])
        );
        assert!(tail.is_empty());
    }

    #[test]
    fn partition_safe_and_tail_retains_last_window_when_min_distance_crosses_chunk() {
        let windows = vec![
            win("chr1", 0, 5, "g1", None),
            win("chr1", 10, 15, "g1", None),
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            Some(4),
            MergeScope::Within,
            Some(0),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        // The first window cannot reach the boundary. Only the suffix remains in tail
        assert_eq!(snapshot(&safe), snapshot(&[win("chr1", 0, 5, "g1", None)]));
        assert_eq!(
            snapshot(&tail),
            snapshot(&[win("chr1", 10, 15, "g1", None)])
        );
    }

    #[test]
    fn partition_safe_and_tail_across_scope_keeps_boundary_suffix_only() {
        let windows = vec![
            win("chr1", 0, 4, "g1", None),
            win("chr1", 5, 7, "g2", None),
            win("chr1", 20, 25, "g1", None),
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            Some(3),
            MergeScope::Across,
            Some(2),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        assert_eq!(safe.len(), 2);
        assert_eq!(snapshot(&tail).len(), 1);
    }

    #[test]
    fn partition_safe_and_tail_across_scope_keeps_overlap_chain_in_tail() {
        let windows = vec![
            // Early windows that should finalize
            win("chr1", 0, 4, "g1", None),
            win("chr1", 15, 18, "g2", None),
            // Overlap/merge chain near the boundary - each link extends reach
            win("chr1", 40, 45, "g3", None),
            win("chr1", 47, 49, "g4", None), // within margin of previous end
            win("chr1", 52, 54, "g5", None), // extends chain to boundary so all three stay in tail
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::Across,
            Some(3),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        assert_eq!(snapshot(&safe).len(), 2);
        assert_eq!(snapshot(&tail).len(), 3);
    }
}

mod tests_mergers {
    use cfdnalab::commands::prepare_windows::{
        config::{CoordinateSet, MergeLabel, MergeScope},
        labels::LabelTuple,
        mergers::{merge_across_groups, merge_windows, merge_within_groups},
        prepare_windows::Window,
    };
    use std::sync::Arc;

    fn win(chrom: &str, start: u32, end: u32, group: &str) -> Window {
        Window {
            chrom: Arc::<str>::from(chrom.to_string()),
            original_start: start,
            original_end: end,
            resized_start: start,
            resized_end: end,
            label_tuples: vec![LabelTuple::new(group.to_string())],
            group_key: group.to_string(),
            score: None,
            merged: false,
        }
    }

    fn snapshot(windows: &[Window]) -> Vec<(String, u32, u32, Vec<String>)> {
        windows
            .iter()
            .map(|w| {
                (
                    w.chrom.as_ref().to_string(),
                    w.resized_start,
                    w.resized_end,
                    w.label_tuples.iter().map(|t| t.input.clone()).collect(),
                )
            })
            .collect()
    }

    #[test]
    fn merge_within_groups_merges_overlaps() {
        let windows = vec![
            win("chr1", 0, 5, "A"),
            win("chr1", 4, 8, "A"),
            win("chr1", 20, 25, "B"),
        ];
        let merged =
            merge_within_groups(windows, 2, MergeLabel::Join, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![
                ("chr1".into(), 0, 8, vec!["A".into()]),
                ("chr1".into(), 20, 25, vec!["B".into()]),
            ]
        );
    }

    #[test]
    fn merge_within_groups_respects_gap_threshold() {
        let windows = vec![win("chr1", 0, 4, "A"), win("chr1", 7, 10, "A")];
        let merged = merge_within_groups(
            windows.clone(),
            2,
            MergeLabel::Join,
            CoordinateSet::Resized,
            false,
        );
        assert_eq!(snapshot(&merged), snapshot(&windows));
    }

    #[test]
    fn merge_within_groups_bridges_gap_within_threshold() {
        let windows = vec![win("chr1", 0, 4, "A"), win("chr1", 6, 9, "A")];
        let merged =
            merge_within_groups(windows, 2, MergeLabel::Join, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 0, 9, vec!["A".into()])]
        );
    }

    #[test]
    fn merge_across_groups_joins_labels() {
        let windows = vec![win("chr1", 0, 4, "G1"), win("chr1", 3, 6, "G2")];
        let merged =
            merge_across_groups(windows, 1, MergeLabel::Join, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 0, 6, vec!["G1".into(), "G2".into()])]
        );
    }

    #[test]
    fn merge_across_groups_sorts_unsorted_input() {
        let windows = vec![win("chr1", 5, 7, "B"), win("chr1", 2, 6, "A")];
        let merged =
            merge_across_groups(windows, 1, MergeLabel::First, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 2, 7, vec!["A".into()])]
        );
    }

    #[test]
    fn merge_across_groups_honors_first_label_policy() {
        let windows = vec![win("chr1", 0, 4, "G1"), win("chr1", 3, 6, "G2")];
        let merged =
            merge_across_groups(windows, 1, MergeLabel::First, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 0, 6, vec!["G1".into()])]
        );
    }

    #[test]
    fn merge_windows_respects_scope_none() {
        let windows = vec![win("chr1", 0, 5, "A")];
        let merged = merge_windows(
            windows.clone(),
            MergeScope::None,
            Some(3),
            MergeLabel::Join,
            CoordinateSet::Resized,
            false,
        );
        assert_eq!(snapshot(&merged), snapshot(&windows));
    }

    #[test]
    fn merge_windows_returns_original_when_gap_missing() {
        let windows = vec![win("chr1", 0, 5, "A"), win("chr1", 10, 12, "A")];
        let merged = merge_windows(
            windows.clone(),
            MergeScope::Within,
            None,
            MergeLabel::Join,
            CoordinateSet::Resized,
            false,
        );
        assert_eq!(snapshot(&merged), snapshot(&windows));
    }
}

mod tests_resizers {
    use cfdnalab::commands::prepare_windows::{
        config::{OobPolicy, PrepareConfig},
        resizers::apply_size_transform,
    };

    fn base_config() -> PrepareConfig {
        let mut cfg = PrepareConfig::default();
        cfg.oob = OobPolicy::Allow;
        cfg
    }

    #[test]
    fn resize_with_odd_size_centers_window() {
        let mut cfg = base_config();
        cfg.resize = Some(5);
        // Odd input length and odd target size yield a single centered placement
        let transformed = apply_size_transform(10, 21, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (21 - 10) / 2 = 15, size 5 spans 13-18
        assert_eq!(transformed, Some((13, 18)));
    }

    #[test]
    fn resize_with_even_size_centers_window_when_parity_matches() {
        let mut cfg = base_config();
        cfg.resize = Some(6);
        // Even input length and even target size yield a single centered placement
        let transformed = apply_size_transform(10, 22, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (22 - 10) / 2 = 16, size 6 spans 13-19
        assert_eq!(transformed, Some((13, 19)));
    }

    #[test]
    fn resize_with_even_input_and_odd_target_picks_left_or_right() {
        let mut cfg = base_config();
        cfg.resize = Some(3);
        // Even input length and odd target size have two equally centered placements
        let transformed = apply_size_transform(10, 16, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (16 - 10) / 2 = 13, size 3 can be 11-14 or 12-15
        let left = (11, 14);
        let right = (12, 15);
        assert!(transformed == Some(left) || transformed == Some(right));

        // TODO: Add examples (different seeds) that shows each outcome (regression tests)
    }

    #[test]
    fn resize_with_odd_input_and_even_target_picks_left_or_right() {
        let mut cfg = base_config();
        cfg.resize = Some(4);
        // Odd input length and even target size have two equally centered placements
        let transformed = apply_size_transform(10, 15, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (15 - 10) / 2 = 12, size 4 can be 10-14 or 11-15
        let left = (10, 14);
        let right = (11, 15);
        assert!(transformed == Some(left) || transformed == Some(right));

        // TODO: Add examples (different seeds) that shows each outcome (regression tests)
    }

    #[test]
    fn flank_with_trim_clamps_to_chrom_bounds() {
        let mut cfg = base_config();
        cfg.flank = Some(vec![5, 5]);
        cfg.oob = OobPolicy::Trim;
        let transformed = apply_size_transform(3, 5, Some(10), &cfg).expect("trim");
        assert_eq!(transformed, Some((0, 10)));
    }

    #[test]
    fn flank_with_drop_returns_none_when_out_of_bounds() {
        let mut cfg = base_config();
        cfg.flank = Some(vec![5, 5]);
        cfg.oob = OobPolicy::Drop;
        let transformed = apply_size_transform(3, 5, Some(6), &cfg).expect("drop");
        assert!(transformed.is_none());
    }

    #[test]
    fn flank_allow_drops_underflow() {
        let mut cfg = base_config();
        cfg.flank = Some(vec![10, 0]);
        cfg.oob = OobPolicy::Allow;
        let transformed = apply_size_transform(2, 4, Some(50), &cfg).expect("allow");
        assert!(transformed.is_none());
    }

    #[test]
    fn trim_policy_returns_none_when_interval_collapses() {
        let mut cfg = base_config();
        cfg.oob = OobPolicy::Trim;
        let transformed = apply_size_transform(10, 11, Some(10), &cfg).expect("trim collapse");
        assert!(transformed.is_none());
    }

    #[test]
    fn no_transform_returns_original_without_bounds_checks() {
        // Arrange
        let mut cfg = base_config();
        cfg.oob = OobPolicy::Drop;

        // Act
        let transformed = apply_size_transform(5, 15, None, &cfg).expect("no transform");

        // Assert
        assert_eq!(transformed, Some((5, 15)));
    }
}

mod tests_parsers {
    use anyhow::Result;
    use cfdnalab::commands::prepare_windows::parsers::{
        parse_cols_indices, parse_distance_bins, parse_record_line, parse_score_filter,
        resolve_column_indices,
    };

    #[test]
    fn parse_distance_bins_and_match_labels() -> Result<()> {
        let bins = parse_distance_bins(&[
            "prox:<10".to_string(),
            "mid:10-20".to_string(),
            "dist:>20".to_string(),
        ])?;
        assert_eq!(bins.match_label(5), Some("prox"));
        assert_eq!(bins.match_label(15), Some("mid"));
        assert_eq!(bins.match_label(50), Some("dist"));
        assert_eq!(bins.match_label(-5), Some("prox"));
        Ok(())
    }

    #[test]
    fn parse_distance_bins_errors_on_invalid_expr() {
        let err = parse_distance_bins(&["bad".to_string()]).unwrap_err();
        assert!(format!("{err}").contains("Invalid distance bin spec"));
    }

    #[test]
    fn parse_distance_bins_prefers_first_matching_rule() -> Result<()> {
        let bins = parse_distance_bins(&["first:<=10".to_string(), "second:<=20".to_string()])?;
        assert_eq!(bins.match_label(5), Some("first"));
        assert_eq!(bins.match_label(15), Some("second"));
        Ok(())
    }

    #[test]
    fn parse_score_filter_evaluates_condition() -> Result<()> {
        let filter = parse_score_filter(">=1.5")?;
        assert!(filter.eval(2.0));
        assert!(!filter.eval(1.0));
        Ok(())
    }

    #[test]
    fn parse_score_filter_errors_on_invalid_operator() {
        let err = parse_score_filter("~=1.0").unwrap_err();
        assert!(format!("{err}").contains("Invalid score filter"));
    }

    #[test]
    fn resolve_indices_and_parse_record_line() -> Result<()> {
        let cols = resolve_column_indices("chrom=0,start=1,end=2", &["3".to_string()], Some("4"))?;
        let (chrom, start, end, group, score) =
            parse_record_line("chr1\t5\t10\tG\t3.5", '\t', &cols)?;
        assert_eq!(chrom, "chr1");
        assert_eq!((start, end), (5, 10));
        assert_eq!(group, "G");
        assert_eq!(score, Some(3.5));
        Ok(())
    }

    #[test]
    fn parse_record_line_handles_missing_group_columns() -> Result<()> {
        let cols = resolve_column_indices("chrom=0,start=1,end=2", &[], None)?;
        let (chrom, start, end, group, score) = parse_record_line("chr1\t0\t5", '\t', &cols)?;
        assert_eq!(chrom, "chr1");
        assert_eq!((start, end), (0, 5));
        assert!(group.is_empty());
        assert!(score.is_none());
        Ok(())
    }

    #[test]
    fn parse_record_line_errors_on_invalid_interval() {
        let cols = resolve_column_indices("chrom=0,start=1,end=2", &[], None).unwrap();
        let err = parse_record_line("chr1\t10\t5", '\t', &cols).unwrap_err();
        assert!(format!("{err}").contains("End must be greater than start"));
    }

    #[test]
    fn parse_cols_indices_requires_all_fields() {
        let err = parse_cols_indices("chrom=0,start=1").unwrap_err();
        assert!(format!("{err}").contains("cols: missing end="));
    }
}

mod tests_near_file {
    use std::{fs::File, io::Write};

    use cfdnalab::commands::prepare_windows::{
        config::{NearDirection, NearEdge},
        near_file::{
            NearChrom, NearDuplicatesPolicy, NearHit, NearInterval, NearSide, NearestResult,
            Strand, load_near_index, nearest_edge_distance,
        },
    };
    use tempfile::TempDir;

    // TODO: Test strand: Strand::Minus and strand: Strand::Unknown!

    #[test]
    fn load_near_index_parses_groups_no_strand() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let path = dir.path().join("near.tsv");
        let mut file = File::create(&path)?;
        writeln!(file, "chr1\t0\t10\tSiteA")?;
        writeln!(file, "chr1\t20\t30\tSiteB")?;
        let index = load_near_index(
            &path,
            '\t',
            false,
            None,
            Some(&[3]),
            false,
            NearDuplicatesPolicy::Error,
        )?;
        let chrom = index.per_chrom.get("chr1").expect("chrom");
        assert_eq!(chrom.intervals.len(), 2);
        assert_eq!(index.group_id_to_name.len(), 2);
        Ok(())
    }

    #[test]
    fn load_near_index_errors_on_overlap() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let path = dir.path().join("near.tsv");
        let mut file = File::create(&path)?;
        writeln!(file, "chr1\t0\t10")?;
        writeln!(file, "chr1\t5\t12")?;
        let err = load_near_index(
            &path,
            '\t',
            false,
            None,
            None,
            false,
            NearDuplicatesPolicy::Error,
        )
        .unwrap_err();
        assert!(format!("{err}").contains("intervals overlap"));
        Ok(())
    }

    #[test]
    fn nearest_edge_distance_handles_overlap_and_sign() {
        let mut chrom = NearChrom {
            intervals: vec![NearInterval {
                start: 10,
                end: 20,
                group_id: Some(0),
                strand: Strand::Plus,
            }],
            cursor: 0,
        };
        let overlap = nearest_edge_distance(
            12,
            18,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            true,
        )
        .unwrap();
        assert_eq!(
            overlap,
            NearestResult::Single(NearHit {
                distance: 0,
                group_id: Some(0),
                side: NearSide::Overlap,
            })
        );

        let window_before = nearest_edge_distance(
            0,
            5,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            true,
        )
        .unwrap();
        assert_eq!(
            window_before,
            NearestResult::Single(NearHit {
                distance: 5,
                group_id: Some(0),
                side: NearSide::Downstream,
            })
        );

        let window_after = nearest_edge_distance(
            25,
            30,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            false,
        )
        .unwrap();
        assert_eq!(
            window_after,
            NearestResult::Single(NearHit {
                distance: 5,
                group_id: Some(0),
                side: NearSide::Upstream,
            })
        );
    }

    #[test]
    fn nearest_edge_distance_zero_on_interval_boundary() {
        let mut chrom = NearChrom {
            intervals: vec![NearInterval {
                start: 10,
                end: 20,
                group_id: Some(0),
                strand: Strand::Plus,
            }],
            cursor: 0,
        };
        let on_boundary = nearest_edge_distance(
            20,
            25,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            false,
        )
        .unwrap();
        assert_eq!(
            on_boundary,
            NearestResult::Single(NearHit {
                distance: 0,
                group_id: Some(0),
                side: NearSide::Overlap,
            })
        );
    }

    #[test]
    fn load_near_index_skips_header() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let path = dir.path().join("near.tsv");
        let mut file = File::create(&path)?;
        writeln!(file, "chrom\tstart\tend")?;
        writeln!(file, "chr1\t0\t5")?;
        let index = load_near_index(
            &path,
            '\t',
            true,
            None,
            None,
            false,
            NearDuplicatesPolicy::Error,
        )?;
        assert_eq!(index.per_chrom.get("chr1").unwrap().intervals.len(), 1);
        Ok(())
    }

    #[test]
    fn nearest_edge_distance_returns_none_without_intervals() {
        let mut chrom = NearChrom {
            intervals: vec![],
            cursor: 0,
        };
        assert!(
            nearest_edge_distance(
                0,
                5,
                &mut chrom,
                &NearEdge::Nearest,
                &NearDirection::Both,
                false
            )
            .is_none()
        );
    }

    #[test]
    fn nearest_edge_distance_respects_left_edge_mode() {
        let mut chrom = NearChrom {
            intervals: vec![NearInterval {
                start: 10,
                end: 20,
                group_id: Some(1),
                strand: Strand::Plus,
            }],
            cursor: 0,
        };
        let dist = nearest_edge_distance(
            30,
            35,
            &mut chrom,
            &NearEdge::Left,
            &NearDirection::Both,
            false,
        )
        .unwrap();
        assert_eq!(
            dist,
            NearestResult::Single(NearHit {
                distance: 20,
                group_id: Some(1),
                side: NearSide::Upstream,
            })
        );
    }

    #[test]
    fn nearest_edge_distance_reports_ties_with_sides() {
        let mut chrom = NearChrom {
            intervals: vec![
                NearInterval {
                    start: 0,
                    end: 5,
                    group_id: Some(1),
                    strand: Strand::Plus,
                },
                NearInterval {
                    start: 25,
                    end: 30,
                    group_id: Some(2),
                    strand: Strand::Plus,
                },
            ],
            cursor: 0,
        };

        let result = nearest_edge_distance(
            10,
            20,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            true,
        )
        .unwrap();
        match result {
            NearestResult::Tie(tie) => {
                assert_eq!(
                    tie.upstream,
                    Some(NearHit {
                        distance: -5,
                        group_id: Some(1),
                        side: NearSide::Upstream,
                    })
                );
                assert_eq!(
                    tie.downstream,
                    Some(NearHit {
                        distance: 5,
                        group_id: Some(2),
                        side: NearSide::Downstream,
                    })
                );
            }
            other => panic!("expected tie, got {other:?}"),
        }
    }
}

mod tests_writers {
    use std::{fs, io::Read};

    use cfdnalab::commands::prepare_windows::{
        config::PrepareConfig,
        labels::{AtomicLabelPart, LabelKey, LabelSchema, LabelTuple},
        prepare_windows::Window,
        writers::{
            ChromTempWriter, concatenate_temps, ensure_temp_writer_for_chrom,
            finalize_temp_writers, write_windows,
        },
    };
    use fxhash::FxHashMap;
    use tempfile::TempDir;

    fn win(chrom: &str, start: u32, end: u32, group: &str) -> Window {
        Window {
            chrom: chrom.to_string().into(),
            original_start: start,
            original_end: end,
            resized_start: start,
            resized_end: end,
            label_tuples: vec![LabelTuple::new(group.to_string())],
            group_key: group.to_string(),
            score: None,
            merged: false,
        }
    }

    fn label_schema() -> LabelSchema {
        LabelSchema::new(&[]).expect("label schema")
    }

    fn out_labels() -> Vec<LabelKey> {
        vec![LabelKey::Atomic(AtomicLabelPart::Input)]
    }

    #[test]
    fn write_windows_outputs_expected_columns() {
        let windows = vec![win("chr1", 0, 5, ""), win("chr1", 10, 15, "grp")];
        let mut buf = Vec::new();
        let schema = label_schema();
        let labels = out_labels();
        write_windows(&mut buf, &windows, '\t', &labels, &schema).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "chr1\t0\t5\t\nchr1\t10\t15\tgrp\n");
    }

    #[test]
    fn write_windows_honors_custom_sep() {
        let windows = vec![win("chr1", 0, 5, "")];
        let mut buf = Vec::new();
        let schema = label_schema();
        let labels = out_labels();
        write_windows(&mut buf, &windows, ',', &labels, &schema).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "chr1,0,5,\n");
    }

    #[test]
    fn ensure_temp_writer_creates_and_reuses_writer() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
        {
            let writer = ensure_temp_writer_for_chrom("chr/1", dir.path(), &mut writers)?;
            let schema = label_schema();
            let labels = out_labels();
            write_windows(
                writer.writer(),
                &[win("chr1", 0, 5, "")],
                '\t',
                &labels,
                &schema,
            )?;
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
        assert_eq!(filename, "chrom.chr_1.bed.tmp");
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
            let schema = label_schema();
            let labels = out_labels();
            write_windows(
                writer.writer(),
                &[win("chr1", 0, 5, "")],
                '\t',
                &labels,
                &schema,
            )?;
        }
        let entries = finalize_temp_writers(&mut writers)?;
        assert!(writers.is_empty());
        assert_eq!(entries.len(), 1);
        assert!(entries[0].1.exists());
        Ok(())
    }

    #[test]
    fn concatenate_temps_writes_all_groups() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let temp_path = dir.path().join("chr1.tmp");
        fs::write(
            &temp_path,
            "chr1\t0\t5\tG\nchr1\t10\t15\tG\nchr1\t20\t25\tH\n",
        )?;
        let mut cfg = PrepareConfig::default();
        cfg.output = dir.path().join("out.tsv");
        cfg.sep = '\t';
        cfg.group_cols = vec!["3".to_string()];
        concatenate_temps(&cfg, &[("chr1".to_string(), temp_path)])?;
        let mut output = String::new();
        fs::File::open(cfg.output)?.read_to_string(&mut output)?;
        assert!(output.contains("chr1\t0\t5\tG"));
        assert!(output.contains("chr1\t20\t25\tH"));
        Ok(())
    }

    #[test]
    fn concatenate_temps_without_groups_writes_three_columns() -> anyhow::Result<()> {
        let dir = TempDir::new()?;
        let temp_path = dir.path().join("chr1.tmp");
        fs::write(&temp_path, "chr1\t0\t5\nchr1\t10\t15\n")?;
        let mut cfg = PrepareConfig::default();
        cfg.output = dir.path().join("out.tsv");
        cfg.sep = '\t';
        concatenate_temps(&cfg, &[("chr1".to_string(), temp_path)])?;
        let mut output = String::new();
        fs::File::open(cfg.output)?.read_to_string(&mut output)?;
        assert_eq!(output.trim(), "chr1\t0\t5\nchr1\t10\t15".trim());
        Ok(())
    }
}

mod tests_chunk {
    use cfdnalab::commands::prepare_windows::{
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
        Window {
            chrom: Arc::<str>::from(chrom.to_string()),
            original_start: start,
            original_end: end,
            resized_start: start,
            resized_end: end,
            label_tuples: vec![LabelTuple::new(group.to_string())],
            group_key: group.to_string(),
            score: None,
            merged: false,
        }
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
        assert_eq!(parsed.start, 0);
        assert_eq!(parsed.end, 10);
        assert_eq!(inputs, vec!["g1", "g2"]);
        Ok(())
    }
}

#[cfg(unix)]
mod tests_stdio {
    use anyhow::{Result, anyhow};
    use cfdnalab::commands::prepare_windows::{
        config::{HeaderMode, OobPolicy, PrepareConfig},
        prepare_windows::run,
    };
    use libc;
    use std::{
        io::Read,
        os::fd::{FromRawFd, RawFd},
        path::PathBuf,
    };

    fn run_with_piped_stdio<F>(input: &str, f: F) -> Result<String>
    where
        F: FnOnce() -> Result<()>,
    {
        unsafe {
            let mut in_pipe: [RawFd; 2] = [0; 2];
            let mut out_pipe: [RawFd; 2] = [0; 2];
            if libc::pipe(in_pipe.as_mut_ptr()) == -1 {
                return Err(anyhow!("pipe failed for stdin"));
            }
            if libc::pipe(out_pipe.as_mut_ptr()) == -1 {
                return Err(anyhow!("pipe failed for stdout"));
            }

            // preload stdin pipe with input
            let bytes = input.as_bytes();
            let mut written = 0usize;
            while written < bytes.len() {
                let chunk = libc::write(
                    in_pipe[1],
                    bytes[written..].as_ptr() as *const _,
                    (bytes.len() - written) as _,
                );
                if chunk <= 0 {
                    return Err(anyhow!("write to stdin pipe failed"));
                }
                written += chunk as usize;
            }
            libc::close(in_pipe[1]);

            let stdin_backup = libc::dup(libc::STDIN_FILENO);
            if stdin_backup == -1 {
                return Err(anyhow!("dup stdin failed"));
            }
            let stdout_backup = libc::dup(libc::STDOUT_FILENO);
            if stdout_backup == -1 {
                return Err(anyhow!("dup stdout failed"));
            }

            if libc::dup2(in_pipe[0], libc::STDIN_FILENO) == -1 {
                return Err(anyhow!("dup2 stdin failed"));
            }
            libc::close(in_pipe[0]);

            if libc::dup2(out_pipe[1], libc::STDOUT_FILENO) == -1 {
                return Err(anyhow!("dup2 stdout failed"));
            }
            libc::close(out_pipe[1]);

            let result = f();
            libc::fflush(std::ptr::null_mut());

            if libc::dup2(stdout_backup, libc::STDOUT_FILENO) == -1 {
                return Err(anyhow!("restore stdout failed"));
            }
            libc::close(stdout_backup);

            if libc::dup2(stdin_backup, libc::STDIN_FILENO) == -1 {
                return Err(anyhow!("restore stdin failed"));
            }
            libc::close(stdin_backup);

            let mut output_file = std::fs::File::from_raw_fd(out_pipe[0]);
            let mut output = String::new();
            output_file.read_to_string(&mut output)?;

            result?;
            Ok(output)
        }
    }

    #[test]
    fn run_supports_stdio() -> Result<()> {
        let mut cfg = PrepareConfig::default();
        cfg.input = PathBuf::from("-");
        cfg.output = PathBuf::from("-");
        cfg.header = HeaderMode::Absent;
        cfg.oob = OobPolicy::Allow;
        cfg.resize = Some(8);
        let chrom_sizes = NamedTempFile::new()?;
        fs::write(chrom_sizes.path(), "chr1\t1000\n")?;
        cfg.chrom_sizes = Some(chrom_sizes.path().to_path_buf());

        let input = "chr1\t0\t5\nchr1\t5\t10\n";
        let output = run_with_piped_stdio(input, || run(&cfg))?;
        let lines: Vec<&str> = output
            .lines()
            .filter(|line| line.starts_with("chr"))
            .collect();
        assert_eq!(lines, vec!["chr1\t0\t6\t", "chr1\t4\t12\t"]);
        Ok(())
    }
}
