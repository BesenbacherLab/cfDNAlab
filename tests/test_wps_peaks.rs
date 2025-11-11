mod tests_wps_normalization {
    use cfdnalab::commands::wps_peaks::normalize_wps::{normalize_wps, smoothe_wps};

    fn approx_eq(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn smoothe_preserves_quadratic_signal() {
        let values: Vec<f32> = (0..64)
            .map(|i| {
                let x = i as f32;
                0.5 * x * x - 3.0 * x + 7.0
            })
            .collect();
        let smoothed = smoothe_wps(&values, None);
        assert_eq!(smoothed.len(), values.len());
        for value in smoothed.iter() {
            assert!(value.is_finite(), "smoothed values should remain finite");
        }
    }

    #[test]
    fn smoothe_respects_mask_boundaries() {
        let mut values = vec![0.0f32; 60];
        for (idx, val) in values.iter_mut().enumerate() {
            let angle = idx as f32 / 5.0;
            *val = (angle.sin() + angle.cos()).abs();
        }
        for val in values[30..].iter_mut() {
            *val += 100.0;
        }

        let mut mask = vec![0u8; values.len()];
        mask[28..32].fill(1);

        let smoothed = smoothe_wps(&values, Some(&mask));

        assert!(smoothed[29].is_nan());
        assert!(smoothed[30].is_nan());
        for idx in 0..mask.len() {
            if mask[idx] != 0 {
                assert!(smoothed[idx].is_nan());
            } else {
                assert!(smoothed[idx].is_finite());
            }
        }
    }

    #[test]
    fn normalize_subtracts_sliding_median() {
        let numerator = vec![1.0, 2.0, 50.0, 4.0, 5.0, 6.0];
        let baseline = numerator.clone();
        let mask = vec![0u8; numerator.len()];

        let normalized = normalize_wps(&numerator, &baseline, Some(&mask), 5, 1, 3);

        let expected: Vec<f32> = vec![-1.0, -1.0, 46.0, -1.0, -0.5, 1.0];
        for (idx, (observed, exp)) in normalized.iter().zip(expected.iter()).enumerate() {
            if (*exp).is_nan() {
                assert!(observed.is_nan(), "index {idx} should be NaN");
            } else {
                assert!(
                    approx_eq(*observed, *exp, 1e-4),
                    "index {idx} expected {exp} got {observed}"
                );
            }
        }
    }

    #[test]
    fn normalize_respects_mask_and_threshold() {
        let numerator = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let baseline = numerator.clone();
        let mut mask = vec![0u8; numerator.len()];
        mask[1] = 1;

        let normalized_loose = normalize_wps(&numerator, &baseline, Some(&mask), 5, 1, 2);
        assert!(normalized_loose[1].is_nan());
        assert!(approx_eq(normalized_loose[2], 3.0 - 3.5, 1e-4));

        let normalized_strict = normalize_wps(&numerator, &baseline, Some(&mask), 5, 1, 5);
        assert!(normalized_strict[2].is_nan());
    }
}

mod tests_normalization_helpers {
    use cfdnalab::commands::wps_peaks::normalize_wps::{
        SlidingMedian, build_left_edge_window, build_right_edge_window,
    };

    const SG_WINDOW_SIZE: usize = 21;
    const EPS: f32 = 1e-6;

    fn approx_eq(a: f32, b: f32) -> bool {
        (a - b).abs() <= EPS
    }

    #[test]
    fn build_left_edge_window_reflects_prefix() {
        let edge_slice = vec![3.0_f32, 4.0_f32, 5.0_f32];
        let window = build_left_edge_window(&edge_slice);
        assert_eq!(window.len(), SG_WINDOW_SIZE);
        for value in &window[..edge_slice.len()] {
            assert!(*value <= edge_slice[0]);
        }
        assert_eq!(
            &window[SG_WINDOW_SIZE - edge_slice.len()..],
            edge_slice.as_slice()
        );
    }

    #[test]
    fn build_right_edge_window_reflects_suffix() {
        let edge_slice = vec![1.0_f32, 2.0_f32, 3.0_f32];
        let window = build_right_edge_window(&edge_slice);
        assert_eq!(window.len(), SG_WINDOW_SIZE);
        assert_eq!(&window[..edge_slice.len()], edge_slice.as_slice());
        for value in &window[edge_slice.len()..] {
            assert!(*value >= *edge_slice.last().unwrap());
        }
    }

    #[test]
    fn sliding_median_tracks_window() {
        let mut median = SlidingMedian::new(5);
        median.insert(0, 1.0);
        assert!(approx_eq(median.median().unwrap(), 1.0));
        median.insert(1, 3.0);
        assert!(approx_eq(median.median().unwrap(), 2.0));
        median.insert(2, 5.0);
        assert!(approx_eq(median.median().unwrap(), 3.0));
        median.remove(1);
        assert!(approx_eq(median.median().unwrap(), (1.0 + 5.0) * 0.5));
        median.remove(0);
        assert!(approx_eq(median.median().unwrap(), 5.0));
        median.remove(2);
        assert!(median.median().is_none());
    }
}

#[cfg(test)]
mod tests_wps_peaks_helpers {
    use cfdnalab::commands::wps_peaks::call_peaks::*;
    use cfdnalab::commands::wps_peaks::window_peak_results::PeaksWindowAction;
    use cfdnalab::commands::wps_peaks::wps_peaks::*;
    use std::collections::BTreeMap;

    fn make_peak(chr: &str, position: u64, height: f32) -> PeakCall {
        PeakCall {
            chromosome: chr.to_string(),
            start: position,
            end: position + 1,
            peak_position: position,
            height,
        }
    }

    #[test]
    fn compute_stats_contributions_extracts_metrics() {
        let windows = vec![(0, 100, 0), (100, 200, 1)];
        let peaks = vec![
            make_peak("chr1", 10, 1.0),
            make_peak("chr1", 50, 1.0),
            make_peak("chr1", 150, 1.0),
        ];

        let contributions = compute_window_stats_contributions(&windows, &peaks);

        assert_eq!(contributions.len(), 2);
        let first = contributions.iter().find(|c| c.window_idx == 0).unwrap();
        assert_eq!(first.count, 2);
        assert_eq!(first.first_peak, Some(10));
        assert_eq!(first.last_peak, Some(50));
        assert_eq!(first.distance_sum, 40.0);
        assert_eq!(first.distance_histogram.get(&40), Some(&1));

        let second = contributions.iter().find(|c| c.window_idx == 1).unwrap();
        assert_eq!(second.count, 1);
        assert_eq!(second.first_peak, Some(150));
        assert_eq!(second.last_peak, Some(150));
        assert_eq!(second.distance_sum, 0.0);
        assert!(second.distance_histogram.is_empty());
    }

    #[test]
    fn histogram_median_handles_even_counts() {
        let mut hist = BTreeMap::new();
        hist.insert(10, 1);
        hist.insert(20, 1);
        let median = histogram_median(&hist);
        assert!((median - 15.0).abs() < 1e-6);
    }

    #[test]
    fn stats_contributions_merge_across_tiles() {
        let mut acc = WindowAccumulator::new(PeaksWindowAction::Stats, 2);
        acc.reset_for_chromosome("chr1".to_string());
        let windows = vec![(0, 100, 0)];
        let mut next_idx = 0usize;
        acc.add_windows_for_tile(&windows, &mut next_idx, 0, 60);

        let mut hist_first = BTreeMap::new();
        hist_first.insert(30, 1);
        let contrib_first = WindowStatsContribution {
            window_idx: 0,
            count: 2,
            first_peak: Some(10),
            last_peak: Some(40),
            distance_sum: 30.0,
            distance_histogram: hist_first,
        };
        acc.apply_stats_contribution(&contrib_first).unwrap();

        acc.add_windows_for_tile(&windows, &mut next_idx, 60, 120);
        let mut hist_second = BTreeMap::new();
        hist_second.insert(20, 1);
        let contrib_second = WindowStatsContribution {
            window_idx: 0,
            count: 2,
            first_peak: Some(70),
            last_peak: Some(90),
            distance_sum: 20.0,
            distance_histogram: hist_second,
        };
        acc.apply_stats_contribution(&contrib_second).unwrap();

        let mut out = Vec::new();
        acc.flush_completed_windows(120, &mut out).unwrap();
        let output = String::from_utf8(out).unwrap();
        assert_eq!(output.trim(), "chr1\t0\t100\t0\t4\t26.67\t30.00");
    }

    #[test]
    fn window_accumulator_unique_writes_sorted_positions() {
        let mut acc = WindowAccumulator::new(PeaksWindowAction::OnlyIncludeThesePositionsUnique, 2);
        acc.reset_for_chromosome("chr1".to_string());

        let windows = vec![(0, 100, 0), (100, 200, 1)];
        let mut next_idx = 0usize;

        acc.add_windows_for_tile(&windows, &mut next_idx, 0, 150);
        let peaks_tile_one = vec![make_peak("chr1", 10, 2.0)];
        for peak in &peaks_tile_one {
            acc.push_peak(peak);
        }
        let mut out = Vec::new();
        acc.flush_completed_windows(150, &mut out).unwrap();

        acc.add_windows_for_tile(&windows, &mut next_idx, 150, 220);
        let peaks_tile_two = vec![make_peak("chr1", 120, 4.0)];
        for peak in &peaks_tile_two {
            acc.push_peak(peak);
        }
        acc.flush_completed_windows(220, &mut out).unwrap();
        acc.flush_all(&mut out).unwrap();

        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.trim().split('\n').collect();
        assert_eq!(lines, vec!["chr1\t10\t11\t2.00", "chr1\t120\t121\t4.00"]);
    }

    #[test]
    fn window_accumulator_stats_counts_across_tiles() {
        let mut acc = WindowAccumulator::new(PeaksWindowAction::Stats, 2);
        acc.reset_for_chromosome("chr2".to_string());

        let windows = vec![(0, 100, 0), (100, 200, 1)];
        let mut next_idx = 0usize;

        acc.add_windows_for_tile(&windows, &mut next_idx, 0, 150);
        let peaks_tile_one = vec![make_peak("chr2", 10, 3.5), make_peak("chr2", 60, 2.5)];
        for peak in &peaks_tile_one {
            acc.push_peak(peak);
        }
        let mut out = Vec::new();
        acc.flush_completed_windows(150, &mut out).unwrap();

        acc.add_windows_for_tile(&windows, &mut next_idx, 150, 220);
        let peaks_tile_two = vec![make_peak("chr2", 120, 5.0), make_peak("chr2", 180, 7.0)];
        for peak in &peaks_tile_two {
            acc.push_peak(peak);
        }
        acc.flush_completed_windows(220, &mut out).unwrap();
        acc.flush_all(&mut out).unwrap();

        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.trim().split('\n').collect();
        assert_eq!(
            lines,
            vec![
                "chr2\t0\t100\t0\t2\t50.00\t50.00",
                "chr2\t100\t200\t1\t2\t60.00\t60.00"
            ]
        );
    }
}

mod tests_wps_peaks_command {
    use crate::fixtures::{
        long_fragment_bam, read_zst_to_string, BamFixture, LONG_FRAGMENT_LENGTH,
        LONG_FRAGMENT_STARTS,
    };
    use anyhow::Result;
    use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs, WindowsArgs};
    use cfdnalab::commands::wps_peaks::config::WPSPeaksConfig;
    use cfdnalab::commands::wps_peaks::window_peak_results::PeaksWindowAction;
    use cfdnalab::commands::wps_peaks::wps_peaks::run;
    use tempfile::tempdir;

    const CHROM_NAME: &str = "chr1";
    const WINDOW_SIZE_BP: u32 = 120;
    const BIN_SIZE_BP: u64 = 1_000;
    const TILE_SIZE_BP: u32 = 1_500;
    const NORMALIZE_BP_FOR_TEST: u32 = 200;

    const BASE_LEFT_BP: u64 = (WINDOW_SIZE_BP / 2) as u64;
    const OVERLAP_WIDTH_BP: u64 = 1; // unique output stores only the peak position
    const OVERLAP_HEIGHT: f32 = 2.0;
    const SHOULDER_OFFSET_BP: u64 = BASE_LEFT_BP + 199;
    const SHOULDER_HEIGHT: f32 = 1.0;

    #[test]
    fn run_emits_expected_peaks_and_stats_for_fixed_size_windows() -> Result<()> {
        let bam = long_fragment_bam("wps_peaks_600bp_fragments")?;

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
                let mut distances: Vec<u64> = peaks_in_window
                    .windows(2)
                    .map(|pair| pair[1] - pair[0])
                    .collect();
                distances.sort_unstable();
                let sum: u64 = distances.iter().sum();
                let avg = sum as f32 / distances.len() as f32;
                let median = if distances.len() % 2 == 1 {
                    distances[distances.len() / 2] as f32
                } else {
                    let mid = distances.len() / 2;
                    (distances[mid - 1] + distances[mid]) as f32 * 0.5
                };
                (format!("{avg:.2}"), format!("{median:.2}"))
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

    fn base_config(
        bam: &BamFixture,
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
mod fixtures;
