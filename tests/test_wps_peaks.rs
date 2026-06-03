#![cfg(all(feature = "cmd_wps_peaks", feature = "testing"))]

mod tests_wps_peaks_command {
    use anyhow::Result;
    use cfdnalab::RunOptions;
    use cfdnalab::run_like_cli::common::{ChromosomeArgs, IOCArgs, WindowsArgs};
    use cfdnalab::run_like_cli::wps_peaks::{PeaksWindowAction, WPSPeaksConfig, run_wps_peaks};
    use cfdnalab::testing::{
        Bed4Row, TempBam, bam_from_fragments, long_inward_fragment_series_bam, read_zst_to_string,
        write_bed4,
    };
    use std::fs::File;
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use tempfile::tempdir;

    const CHROM_NAME: &str = "chr1";
    const LONG_FRAGMENT_LENGTH: i64 = 600;
    const LONG_FRAGMENT_STARTS: [i64; 10] =
        [0, 400, 800, 1_200, 1_600, 2_000, 2_400, 2_800, 3_200, 3_600];
    const WINDOW_SIZE_BP: u32 = 120;
    const BIN_SIZE_BP: u64 = 1_000;
    const TILE_SIZE_BP: u32 = 1_500;
    const NORMALIZE_BP_FOR_TEST: u32 = 200;

    const BASE_LEFT_BP: u64 = (WINDOW_SIZE_BP / 2) as u64;
    const OVERLAP_WIDTH_BP: u64 = 1; // unique output stores only the peak position
    const OVERLAP_HEIGHT: f32 = 2.0;
    const SHOULDER_OFFSET_BP: u64 = BASE_LEFT_BP + 199;
    const SHOULDER_HEIGHT: f32 = 1.0;

    fn run(cfg: &WPSPeaksConfig) -> Result<()> {
        run_wps_peaks(cfg, RunOptions::new_quiet()).map(|_| ())
    }

    fn empty_three_chrom_bam(name: &str) -> Result<TempBam> {
        bam_from_fragments(
            name,
            vec![
                ("chr1".to_string(), 1_000),
                ("chr2".to_string(), 1_000),
                ("chr3".to_string(), 1_000),
            ],
            Vec::new(),
            Vec::new(),
        )
    }

    // KEEP-IN-TESTS: wps-peaks command output or artifact behavior.
    #[test]
    fn run_emits_expected_peaks_and_stats_for_fixed_size_windows() -> Result<()> {
        let bam = long_inward_fragment_series_bam("wps_peaks_600bp_fragments")?;

        // Per-window unique output captures the expected peak coordinates.
        let peaks_dir = tempdir()?;
        let peaks_cfg = base_config(
            &bam,
            peaks_dir.path(),
            "by_size_peaks",
            Some(PeaksWindowAction::OnlyIncludeThesePositionsUnique),
        );
        run(&peaks_cfg)?;
        let peaks_path = peaks_dir
            .path()
            .join("by_size_peaks.wps.peaks.unique.tsv.zst");
        let peak_rows = parse_unique_peaks(&read_zst_to_string(&peaks_path)?);
        // Peaks Derivation:
        // - Each fragment produces a single-fragment shoulder roughly `start + BASE_LEFT_BP + 199`
        //   bases into the insert (where the residual plateau reaches its maximum after the median
        //   subtraction). These shoulders sit at height ~1.0.
        // - Consecutive fragments overlap on `[start_i + BASE_LEFT_BP + 400, start_i + BASE_LEFT_BP + 480)`
        //   so their residuals rise to ~2.0 there. The unique output collapses each Snyder peak to the
        //   position of the maximum, i.e., the left edge of those overlap bands.
        let mut expected_peaks: Vec<PeakRow> = LONG_FRAGMENT_STARTS
            .windows(2)
            .map(|pair| {
                let start = (pair[1] as u64) + BASE_LEFT_BP;
                PeakRow {
                    start,
                    end: start + OVERLAP_WIDTH_BP,
                    height: OVERLAP_HEIGHT,
                }
            })
            .collect();
        expected_peaks.extend(
            LONG_FRAGMENT_STARTS
                .iter()
                .copied()
                .enumerate()
                .skip(1) // First fragment lacks a leading shoulder
                .take(LONG_FRAGMENT_STARTS.len().saturating_sub(2)) // Last fragment has no trailing shoulder
                .map(|(_, start)| {
                    let pos = (start as u64) + SHOULDER_OFFSET_BP;
                    PeakRow {
                        start: pos,
                        end: pos + 1,
                        height: SHOULDER_HEIGHT,
                    }
                }),
        );
        expected_peaks.sort_by_key(|peak| peak.start);
        assert_eq!(peak_rows.len(), expected_peaks.len());
        for (actual, expected) in peak_rows.iter().zip(expected_peaks.iter()) {
            assert_eq!(actual.start, expected.start);
            assert_eq!(actual.end, expected.end);
            assert!(
                (actual.height - expected.height).abs() < 1e-6,
                "expected height {} got {}",
                expected.height,
                actual.height
            );
        }

        // Stats output with fixed-size windows
        let stats_dir = tempdir()?;
        let stats_cfg = base_config(
            &bam,
            stats_dir.path(),
            "by_size_stats",
            Some(PeaksWindowAction::Stats),
        );
        run(&stats_cfg)?;
        let stats_path = stats_dir
            .path()
            .join("by_size_stats.wps.peaks.stats.tsv.zst");
        let mut stats_rows = parse_stats(&read_zst_to_string(&stats_path)?);
        stats_rows.sort_by_key(|row| row.index);
        // Windows: binning the 4.7kb contig yields indices 0-4. Populated windows capture the repeated
        // 400bp spacing between adjacent overlaps, so both average and median distances equal 400bp,
        // while windows with fewer than two peaks report `NaN`.
        let peak_positions: Vec<u64> = expected_peaks.iter().map(|peak| peak.start).collect();
        let chrom_len_bp = LONG_FRAGMENT_STARTS.last().copied().unwrap_or(0) as u64
            + LONG_FRAGMENT_LENGTH as u64
            + 500;
        let window_count = ((chrom_len_bp + BIN_SIZE_BP - 1) / BIN_SIZE_BP) as u64;
        let mut expected_stats = Vec::new();
        for idx in 0..window_count {
            let window_start = idx * BIN_SIZE_BP;
            let window_end = (window_start + BIN_SIZE_BP).min(chrom_len_bp);
            let peaks_in_window: Vec<u64> = peak_positions
                .iter()
                .copied()
                .filter(|pos| *pos >= window_start && *pos < window_end)
                .collect();
            let count = peaks_in_window.len() as u32;
            let (avg_distance, median_distance) = if count < 2 {
                ("NaN".to_string(), "NaN".to_string())
            } else {
                // Stats are aggregated per tile, so windows that straddle the 1.5kb tile boundary only
                // see intra-tile gaps. Drop cross-tile pairs to match the production logic.
                let mut distances: Vec<u64> = peaks_in_window
                    .windows(2)
                    .filter_map(|pair| {
                        let left_tile = pair[0] / TILE_SIZE_BP as u64;
                        let right_tile = pair[1] / TILE_SIZE_BP as u64;
                        (left_tile == right_tile).then_some(pair[1] - pair[0])
                    })
                    .collect();
                distances.sort_unstable();
                if distances.is_empty() {
                    ("NaN".to_string(), "NaN".to_string())
                } else {
                    let sum: u64 = distances.iter().sum();
                    let avg = sum as f32 / distances.len() as f32;
                    let median = if distances.len() % 2 == 1 {
                        distances[distances.len() / 2] as f32
                    } else {
                        let mid = distances.len() / 2;
                        (distances[mid - 1] + distances[mid]) as f32 * 0.5
                    };
                    (format!("{avg:.2}"), format!("{median:.2}"))
                }
            };
            expected_stats.push(StatsRow {
                start: window_start,
                end: window_end,
                index: idx,
                count,
                avg_distance,
                median_distance,
            });
        }
        assert_eq!(stats_rows, expected_stats);

        Ok(())
    }

    // KEEP-IN-TESTS: wps-peaks command output or artifact behavior.
    #[test]
    fn global_mode_handles_three_chromosomes() -> Result<()> {
        let bam = empty_three_chrom_bam("wps_peaks_three_chr_global")?;
        let out_dir = tempdir()?;
        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec![
                "chr1".to_string(),
                "chr2".to_string(),
                "chr3".to_string(),
            ]),
            chromosomes_file: None,
        };
        let mut cfg = WPSPeaksConfig::new(ioc, chromosomes, None);
        cfg.shared_args
            .set_output_prefix("three_chr_global".to_string());
        cfg.shared_args.set_window_size(WINDOW_SIZE_BP);
        cfg.shared_args.set_decimals(2);
        cfg.shared_args.set_tile_size(TILE_SIZE_BP);
        cfg.shared_args.set_min_fragment_length(WINDOW_SIZE_BP);
        cfg.shared_args.set_max_fragment_length(2_000);
        cfg.shared_args.set_min_mapq(0);
        cfg.no_smoothing = true;
        cfg.normalize_bp = NORMALIZE_BP_FOR_TEST;
        cfg.min_unmasked = 10;
        cfg.min_peak_height = 0.75;

        run(&cfg)?;

        let peaks_path = out_dir.path().join("three_chr_global.wps.peaks.tsv.zst");
        let text = read_zst_to_string(&peaks_path)?;
        assert!(
            text.lines()
                .eq(["chromosome\tstart\tend\tpeak_position\theight"]),
            "Empty three-chromosome input should produce only the global peaks header"
        );

        Ok(())
    }

    // KEEP-IN-TESTS: wps-peaks command output or artifact behavior.
    #[test]
    fn by_size_stats_handles_three_chromosomes() -> Result<()> {
        let bam = empty_three_chrom_bam("wps_peaks_three_chr_by_size")?;
        let out_dir = tempdir()?;
        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec![
                "chr1".to_string(),
                "chr2".to_string(),
                "chr3".to_string(),
            ]),
            chromosomes_file: None,
        };
        let mut cfg = WPSPeaksConfig::new(ioc, chromosomes, Some(PeaksWindowAction::Stats));
        cfg.shared_args
            .set_output_prefix("three_chr_by_size".to_string());
        cfg.shared_args.set_window_size(WINDOW_SIZE_BP);
        cfg.shared_args.set_decimals(2);
        cfg.shared_args.set_tile_size(TILE_SIZE_BP);
        cfg.shared_args.set_min_fragment_length(WINDOW_SIZE_BP);
        cfg.shared_args.set_max_fragment_length(2_000);
        cfg.shared_args.set_min_mapq(0);
        cfg.shared_args.set_windows(WindowsArgs {
            by_size: Some(1_000),
            by_bed: None,
        });
        cfg.no_smoothing = true;
        cfg.normalize_bp = NORMALIZE_BP_FOR_TEST;
        cfg.min_unmasked = 10;
        cfg.min_peak_height = 0.75;

        run(&cfg)?;

        let stats_path = out_dir
            .path()
            .join("three_chr_by_size.wps.peaks.stats.tsv.zst");
        let text = read_zst_to_string(&stats_path)?;
        let lines: Vec<_> = text.lines().collect();
        // Fixed-size window indices are global across the requested chromosome order. Each 1kb
        // chromosome has one [0, 1000) bin, so chr1/chr2/chr3 get indices 0/1/2.
        assert_eq!(
            lines,
            vec![
                "chromosome\tstart\tend\twindow_index\tcount\tavg_distance\tmedian_distance",
                "chr1\t0\t1000\t0\t0\tNaN\tNaN",
                "chr2\t0\t1000\t1\t0\tNaN\tNaN",
                "chr3\t0\t1000\t2\t0\tNaN\tNaN",
            ]
        );

        Ok(())
    }

    // KEEP-IN-TESTS: wps-peaks command output or artifact behavior.
    #[test]
    fn by_bed_stats_handles_three_chromosomes() -> Result<()> {
        let bam = empty_three_chrom_bam("wps_peaks_three_chr_by_bed")?;
        let out_dir = tempdir()?;
        let bed_path = out_dir.path().join("three_chr_windows.bed");
        write_bed4(
            &bed_path,
            &[
                Bed4Row::new("chr1", 0, 1000, "chr1_window"),
                Bed4Row::new("chr2", 0, 1000, "chr2_window"),
                Bed4Row::new("chr3", 0, 1000, "chr3_window"),
            ],
        )?;

        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec![
                "chr1".to_string(),
                "chr2".to_string(),
                "chr3".to_string(),
            ]),
            chromosomes_file: None,
        };
        let mut cfg = WPSPeaksConfig::new(ioc, chromosomes, Some(PeaksWindowAction::Stats));
        cfg.shared_args
            .set_output_prefix("three_chr_by_bed".to_string());
        cfg.shared_args.set_window_size(WINDOW_SIZE_BP);
        cfg.shared_args.set_decimals(2);
        cfg.shared_args.set_tile_size(TILE_SIZE_BP);
        cfg.shared_args.set_min_fragment_length(WINDOW_SIZE_BP);
        cfg.shared_args.set_max_fragment_length(2_000);
        cfg.shared_args.set_min_mapq(0);
        cfg.shared_args.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
        });
        cfg.no_smoothing = true;
        cfg.normalize_bp = NORMALIZE_BP_FOR_TEST;
        cfg.min_unmasked = 10;
        cfg.min_peak_height = 0.75;

        run(&cfg)?;

        let stats_path = out_dir
            .path()
            .join("three_chr_by_bed.wps.peaks.stats.tsv.zst");
        let text = read_zst_to_string(&stats_path)?;
        let lines: Vec<_> = text.lines().collect();
        assert_eq!(
            lines,
            vec![
                "chromosome\tstart\tend\twindow_index\tcount\tavg_distance\tmedian_distance",
                "chr1\t0\t1000\t0\t0\tNaN\tNaN",
                "chr2\t0\t1000\t1\t0\tNaN\tNaN",
                "chr3\t0\t1000\t2\t0\tNaN\tNaN",
            ]
        );

        Ok(())
    }

    // KEEP-IN-TESTS: wps-peaks command output or artifact behavior.
    #[test]
    fn blacklist_near_boundary_removes_cross_tile_distance() -> Result<()> {
        let bam = long_inward_fragment_series_bam("wps_peaks_boundary_blacklist")?;
        let bed_dir = tempdir()?;
        let far_blacklist = write_blacklist_file(bed_dir.path(), "far", &[(1_000, 1_200)])?;
        let near_blacklist = write_blacklist_file(bed_dir.path(), "near", &[(2_900, 3_000)])?;

        let far_stats =
            run_stats_with_blacklist(&bam, Some(&far_blacklist), "stats_far", Some(3_000))?;
        let far_mid = far_stats.iter().find(|row| row.index == 3).unwrap();
        assert_ne!(
            far_mid.avg_distance, "NaN",
            "without a boundary mask we should retain the 400bp cross-tile gap"
        );

        let near_stats =
            run_stats_with_blacklist(&bam, Some(&near_blacklist), "stats_near", Some(3_000))?;
        let near_mid = near_stats.iter().find(|row| row.index == 3).unwrap();
        assert_eq!(
            near_mid.avg_distance, "NaN",
            "mask ending at the tile edge should eliminate the cross-tile distance"
        );
        Ok(())
    }

    #[derive(Debug, PartialEq)]
    struct PeakRow {
        start: u64,
        end: u64,
        height: f32,
    }

    #[derive(Debug, PartialEq)]
    struct StatsRow {
        start: u64,
        end: u64,
        index: u64,
        count: u32,
        avg_distance: String,
        median_distance: String,
    }

    fn parse_unique_peaks(text: &str) -> Vec<PeakRow> {
        text.lines()
            .skip(1)
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                PeakRow {
                    start: cols[1].parse().unwrap(),
                    end: cols[2].parse().unwrap(),
                    height: cols[3].parse().unwrap(),
                }
            })
            .collect()
    }

    fn parse_stats(text: &str) -> Vec<StatsRow> {
        text.lines()
            .skip(1)
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                StatsRow {
                    start: cols[1].parse().unwrap(),
                    end: cols[2].parse().unwrap(),
                    index: cols[3].parse().unwrap(),
                    count: cols[4].parse().unwrap(),
                    avg_distance: cols[5].to_string(),
                    median_distance: cols[6].to_string(),
                }
            })
            .collect()
    }

    fn write_blacklist_file(dir: &Path, name: &str, entries: &[(u64, u64)]) -> Result<PathBuf> {
        let path = dir.join(format!("{name}.bed"));
        let mut file = File::create(&path)?;
        for (start, end) in entries {
            writeln!(file, "chr1\t{start}\t{end}")?;
        }
        Ok(path)
    }

    fn run_stats_with_blacklist(
        bam: &TempBam,
        blacklist: Option<&Path>,
        prefix: &str,
        tile_size: Option<u32>,
    ) -> Result<Vec<StatsRow>> {
        let stats_dir = tempdir()?;
        let mut cfg = base_config(
            bam,
            stats_dir.path(),
            prefix,
            Some(PeaksWindowAction::Stats),
        );
        cfg.shared_args.blacklist = blacklist.map(|path| vec![path.to_path_buf()]);
        if let Some(size) = tile_size {
            cfg.set_tile_size(size);
        }
        run(&cfg)?;
        let stats_path = stats_dir
            .path()
            .join(format!("{prefix}.wps.peaks.stats.tsv.zst"));
        let mut stats_rows = parse_stats(&read_zst_to_string(&stats_path)?);
        stats_rows.sort_by_key(|row| row.index);
        Ok(stats_rows)
    }

    fn base_config(
        bam: &TempBam,
        out_dir: &std::path::Path,
        prefix: &str,
        per_window: Option<PeaksWindowAction>,
    ) -> WPSPeaksConfig {
        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.to_path_buf(),
            n_threads: 1,
        };
        let mut chromosomes = ChromosomeArgs::default();
        chromosomes.chromosomes = Some(vec![CHROM_NAME.to_string()]);
        let mut cfg = WPSPeaksConfig::new(ioc, chromosomes, per_window);
        cfg.shared_args.set_output_prefix(prefix.to_string());
        cfg.shared_args.set_window_size(WINDOW_SIZE_BP);
        cfg.shared_args.set_decimals(2);
        cfg.shared_args.set_tile_size(TILE_SIZE_BP);
        cfg.shared_args.set_min_fragment_length(WINDOW_SIZE_BP);
        cfg.shared_args.set_max_fragment_length(2_000);
        cfg.shared_args.set_min_mapq(0);
        cfg.shared_args.set_windows(WindowsArgs {
            by_size: Some(BIN_SIZE_BP),
            ..Default::default()
        });
        cfg.no_smoothing = true;
        cfg.normalize_bp = NORMALIZE_BP_FOR_TEST;
        cfg.min_unmasked = 10;
        // Height 0.75 keeps the overlap-only residuals (height ~1.0) while
        // dropping the single-fragment shoulders (residual ~0.5). This mirrors
        // the manual derivation that only the pairwise overlaps form peaks.
        cfg.min_peak_height = 0.75;
        cfg
    }
}
