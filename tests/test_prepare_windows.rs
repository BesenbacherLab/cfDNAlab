#![cfg(feature = "cmd_prepare_windows")]

mod tests_prepare_windows_pipeline {

    use anyhow::Result;
    use cfdnalab::RunOptions;
    use cfdnalab::run_like_cli::{
        common::BlacklistStrategy,
        prepare_windows::{
            CoordinateSet, DedupKeep, DistSign, DistancePolicy, HeaderMode, MergeLabel, MergeScope,
            NearDirection, NearEdge, NearTiePolicy, OobPolicy, PrepareConfig, run_prepare_windows,
        },
    };
    use flate2::{Compression, write::GzEncoder};
    use std::fs;
    use std::io::{BufWriter, Read, Write};
    use tempfile::TempDir;
    use zstd::Decoder as ZstdDecoder;

    fn run(cfg: &PrepareConfig) -> Result<()> {
        run_prepare_windows(cfg, RunOptions::new_quiet()).map(|_| ())
    }

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
        cfg.blacklist_strategy = BlacklistStrategy::Any;
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
            "win-direction".to_string(),
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
        cfg.blacklist_strategy = BlacklistStrategy::Any;
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
        cfg.out_labels = vec!["win-direction".to_string(), "bin".to_string()];
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;

        assert_eq!(
            lines,
            vec![
                "chr1\t6\t8\t-\tprox".to_string(),
                "chr1\t52\t54\t+\tmid".to_string(),
                "chr1\t150\t152\t+\tfar".to_string(),
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
        cfg.distance_from = CoordinateSet::Resized;
        cfg.flank = Some(vec![5, 5]);
        cfg.out_labels = vec![
            "win-direction".to_string(),
            "near-name".to_string(),
            "bin".to_string(),
        ];
        cfg.oob = OobPolicy::Allow;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        // Resized 15..35 yields a 5 bp distance bin and output keeps resized coordinates
        assert_eq!(lines, vec!["chr1\t15\t35\t+\tSITE\tnear".to_string()]);
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
            "win-direction".to_string(),
            "near-name".to_string(),
            "bin".to_string(),
        ];
        cfg.oob = OobPolicy::Allow;
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        // Original 20..30 yields a 10 bp distance bin while output keeps resized coordinates
        assert_eq!(lines, vec!["chr1\t15\t35\t+\tSITE\tfar".to_string()]);
        Ok(())
    }

    #[test]
    fn should_annotate_ties_with_directional_groups_upstream_downstream() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t10\t20"])?;
        let near = write_temp_file(
            &tmpdir,
            "near.tsv",
            &["chr1\t0\t5\tLEFT", "chr1\t25\t30\tRIGHT"],
        )?;
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
        cfg.out_labels = vec!["win-direction".to_string(), "near-name".to_string()];
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;
        // Downstream side (window right of first interval) is '+'
        // Upstream side (window left of second interval) is '-'
        assert_eq!(lines, vec!["chr1\t10\t20\t+,-\tLEFT,RIGHT".to_string()]);
        Ok(())
    }

    #[test]
    fn should_annotate_ties_with_directional_groups_downstream_downstream() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t10\t20"])?;
        let near = write_temp_file(
            &tmpdir,
            "near.tsv",
            &["chr1\t0\t5\tLEFT\t+", "chr1\t25\t30\tRIGHT\t-"],
        )?;
        let output = tmpdir.path().join("out.tsv");

        let mut cfg = PrepareConfig::default();
        cfg.input = input;
        cfg.output = output;
        cfg.header = HeaderMode::Absent;
        cfg.near = Some(near);
        cfg.near_header = HeaderMode::Absent;
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.near_strand_col = Some("4".to_string());
        cfg.near_edge = NearEdge::Nearest;
        cfg.near_direction = NearDirection::Both;
        cfg.distance_sign = DistSign::Signed;
        cfg.out_labels = vec!["win-direction".to_string(), "near-name".to_string()];
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;
        // Both the left and right near intervals are downstream (different strand-orientations)
        // Ordering follows left of window, then right of window.
        // NOTE: The +,+ is compacted into +
        assert_eq!(lines, vec!["chr1\t10\t20\t+\tLEFT,RIGHT".to_string()]);
        Ok(())
    }

    #[test]
    fn should_drop_ties_when_configured() -> Result<()> {
        let tmpdir = TempDir::new()?;
        let input = write_temp_file(&tmpdir, "input.tsv", &["chr1\t10\t20"])?;
        let near = write_temp_file(
            &tmpdir,
            "near.tsv",
            &["chr1\t0\t5\tLEFT", "chr1\t25\t30\tRIGHT"],
        )?;
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
        cfg.out_labels = vec!["win-direction".to_string(), "near-name".to_string()];
        cfg.oob = OobPolicy::Allow;

        let lines = run_pipeline(&cfg)?;

        assert_eq!(
            lines,
            vec!["chr1\t0\t10\t-\tTARG", "chr1\t150\t160\t+\tTARG"],
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
        cfg.out_labels = vec!["win-direction".to_string()];
        cfg.oob = OobPolicy::Allow;

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        assert_eq!(
            lines,
            vec![
                "chr1\t0\t10\t-".to_string(),
                "chr1\t150\t160\t+".to_string()
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
                "chr2\t10\t22\tD",
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
        cfg.blacklist_strategy = BlacklistStrategy::Any;
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
            "win-direction".to_string(),
            "near-name".to_string(),
            "bin".to_string(),
        ];
        cfg.chrom_sizes = Some(chrom_sizes);

        // Act
        let lines = run_pipeline(&cfg)?;

        // Assert
        // Under OOB allow, the first chr2 window underflows and is dropped. The second survives and is resized.
        // Merge on original coordinates then resize to 10 bp:
        // - chr1 A windows merge to 10..21 then resize to either 10..20 (left) or 11..21 (right) because parity is ambiguous.
        // - chr2 D resizes from 10..22 to 11..21 (unique placement).
        let expected_left = vec![
            "chr1\t10\t20\tA\t=\tN1\tclose".to_string(),
            "chr2\t11\t21\tD\t+\tN2\tclose".to_string(),
        ];
        // When resize parity is ambiguous, placement can shift by one bp depending on the
        // deterministic hash choice (left vs right). Accept either centred placement.
        let expected_right = vec![
            "chr1\t11\t21\tA\t=\tN1\tclose".to_string(),
            "chr2\t11\t21\tD\t+\tN2\tclose".to_string(),
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
