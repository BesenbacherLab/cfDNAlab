mod tests_wps_normalization {
    // MOVE-MODULE-LOCAL: direct private WPS normalization helper tests.

    use crate::commands::wps_peaks::normalize_wps::{normalize_wps, smoothe_wps};

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

        let normalized = normalize_wps(
            &numerator,
            &baseline,
            Some(&mask),
            5, // window_size
            1, // stride
            3, // min_unmasked
        );

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

        let normalized_loose = normalize_wps(
            &numerator,
            &baseline,
            Some(&mask),
            5, // window_size
            1, // stride
            2, // min_unmasked
        );
        assert!(normalized_loose[1].is_nan());
        assert!(approx_eq(normalized_loose[2], 3.0 - 3.5, 1e-4));

        let normalized_strict = normalize_wps(
            &numerator,
            &baseline,
            Some(&mask),
            5, // window_size
            1, // stride
            5, // min_unmasked
        );
        assert!(normalized_strict[2].is_nan());
    }
}

mod tests_normalization_helpers {
    // MOVE-MODULE-LOCAL: direct private WPS normalization helper tests.

    use crate::commands::wps_peaks::normalize_wps::{
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
    // MOVE-MODULE-LOCAL: direct private WPS peak-calling helper tests.

    use crate::commands::wps_peaks::call_peaks::*;
    use crate::commands::wps_peaks::window_peak_results::PeaksWindowAction;
    use crate::commands::wps_peaks::wps_peaks::*;
    use crate::shared::interval::IndexedInterval;
    use crate::shared::tiled_run::Tile;
    use std::collections::BTreeMap;
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use tempfile::NamedTempFile;

    fn make_peak(chr: &str, position: u64, height: f32) -> PeakCall {
        PeakCall::new(chr.to_string(), position, position + 1, position, height, 0)
            .expect("test peak should be valid")
    }

    fn indexed_windows(entries: &[(u64, u64, u64)]) -> Vec<IndexedInterval<u64>> {
        IndexedInterval::from_tuples(entries).expect("test windows should be valid")
    }

    #[test]
    fn compute_stats_contributions_extracts_metrics() {
        let windows = indexed_windows(&[(0, 100, 0), (100, 200, 1)]);
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
    fn stats_distance_summary_returns_nan_without_distances() {
        let hist = BTreeMap::new();
        let (avg, median) = stats_distance_summary(0.0, &hist);
        assert!(avg.is_nan());
        assert!(median.is_nan());
    }

    #[test]
    fn stats_distance_summary_reports_average_and_median() {
        let mut hist = BTreeMap::new();
        hist.insert(25, 2);
        hist.insert(50, 1);
        let (avg, median) = stats_distance_summary(100.0, &hist);
        let expected_avg = (100.0 / 3.0) as f32;
        assert!(
            (avg - expected_avg).abs() < 1e-6,
            "avg {avg} expected {expected_avg}"
        );
        assert_eq!(median, 25.0);
    }

    #[test]
    fn stats_contributions_merge_across_tiles() {
        let mut acc = WindowAccumulator::new(PeaksWindowAction::Stats, 2);
        acc.reset_for_chromosome("chr1".to_string());
        let windows = indexed_windows(&[(0, 100, 0)]);
        let mut next_idx = 0usize;
        acc.add_windows_for_tile(&windows, &mut next_idx, 0, 60);

        let mut hist_first = BTreeMap::new();
        hist_first.insert(30, 1);
        let contrib_first = WindowStatsContribution {
            window_idx: 0,
            count: 2,
            first_peak: Some(10),
            last_peak: Some(40),
            first_segment: Some(0),
            last_segment: Some(0),
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
            first_segment: Some(0),
            last_segment: Some(0),
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

        let windows = indexed_windows(&[(0, 100, 0), (100, 200, 1)]);
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

        let windows = indexed_windows(&[(0, 100, 0), (100, 200, 1)]);
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

    #[test]
    fn stats_contributions_skip_blacklisted_gap_distances() {
        // Window spans 0-200 and collects stats via contributions instead of streaming peaks.
        // First contribution represents a peak before the blacklist, second one after it. Their
        // segment markers differ, so merging must not invent a 100bp gap between them.
        let mut acc = WindowAccumulator::new(PeaksWindowAction::Stats, 2);
        acc.reset_for_chromosome("chr1".to_string());
        let windows = indexed_windows(&[(0, 200, 0)]);
        let mut next_idx = 0usize;
        acc.add_windows_for_tile(&windows, &mut next_idx, 0, 200);

        let contrib_first = WindowStatsContribution {
            window_idx: 0,
            count: 1,
            first_peak: Some(50),
            last_peak: Some(50),
            first_segment: Some(0),
            last_segment: Some(0),
            distance_sum: 0.0,
            distance_histogram: BTreeMap::new(),
        };
        acc.apply_stats_contribution(&contrib_first).unwrap();

        let contrib_second = WindowStatsContribution {
            window_idx: 0,
            count: 1,
            first_peak: Some(150),
            last_peak: Some(150),
            first_segment: Some(85),
            last_segment: Some(85),
            distance_sum: 0.0,
            distance_histogram: BTreeMap::new(),
        };
        acc.apply_stats_contribution(&contrib_second).unwrap();

        let mut out = Vec::new();
        acc.flush_all(&mut out).unwrap();
        let output = String::from_utf8(out).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(fields[5], "NaN");
        assert_eq!(fields[6], "NaN");
    }

    #[test]
    fn aligned_and_buffered_unique_outputs_match() {
        let data = fixed_size_test_data();
        let peak_file = write_peaks_file(&data.peaks);

        let buffered = buffered_unique_rows(
            &data.tile,
            data.windows.as_slice(),
            peak_file.path(),
            data.decimals,
        );
        let aligned = aligned_unique_rows(data.tile.chr.as_str(), peak_file.path(), data.decimals);

        assert_eq!(aligned.trim_end(), buffered.trim_end());
    }

    #[test]
    fn aligned_and_buffered_stats_outputs_match() {
        let data = fixed_size_test_data();
        let peak_file = write_peaks_file(&data.peaks);
        let contributions =
            compute_window_stats_contributions(data.windows.as_slice(), &data.peaks);

        let buffered = buffered_stats_rows(
            &data.tile,
            data.windows.as_slice(),
            peak_file.path(),
            data.decimals,
        );
        let aligned = aligned_stats_rows(
            data.tile.chr.as_str(),
            data.windows.as_slice(),
            contributions.as_slice(),
            data.decimals,
        );

        assert_eq!(aligned.trim_end(), buffered.trim_end());
    }

    #[test]
    fn aligned_and_buffered_unique_outputs_match_across_tiles() {
        let data = two_tile_test_data();
        let peak_files = write_peak_files(&data.peaks_by_tile);
        let paths: Vec<PathBuf> = peak_files
            .iter()
            .map(|file| file.path().to_path_buf())
            .collect();

        let buffered = buffered_unique_rows_multi(
            &data.tiles,
            data.all_windows.as_slice(),
            paths.as_slice(),
            data.decimals,
        );
        let aligned = aligned_unique_rows_multi(&data.tiles, paths.as_slice(), data.decimals);

        assert_eq!(aligned.trim_end(), buffered.trim_end());
    }

    #[test]
    fn aligned_and_buffered_stats_outputs_match_across_tiles() {
        let data = two_tile_test_data();
        let peak_files = write_peak_files(&data.peaks_by_tile);
        let paths: Vec<PathBuf> = peak_files
            .iter()
            .map(|file| file.path().to_path_buf())
            .collect();

        let buffered = buffered_stats_rows_multi(
            &data.tiles,
            data.all_windows.as_slice(),
            paths.as_slice(),
            data.decimals,
        );

        let aligned = aligned_stats_rows_multi(
            &data.tiles,
            data.bin_size,
            data.chrom_len,
            &data.peaks_by_tile,
            data.decimals,
        );

        assert_eq!(aligned.trim_end(), buffered.trim_end());
    }

    struct FixedSizeTestData {
        tile: Tile,
        peaks: Vec<PeakCall>,
        windows: Vec<IndexedInterval<u64>>,
        decimals: usize,
    }

    fn make_tile(
        chr: &str,
        index: u32,
        core_start: u32,
        core_end: u32,
        fetch_start: u32,
        fetch_end: u32,
    ) -> Tile {
        Tile::from_coords(
            chr.to_string(),
            0,
            index,
            core_start,
            core_end,
            fetch_start,
            fetch_end,
        )
        .expect("test tile should be valid")
    }

    struct MultiTileTestData {
        tiles: Vec<Tile>,
        peaks_by_tile: Vec<Vec<PeakCall>>,
        all_windows: Vec<IndexedInterval<u64>>,
        bin_size: u64,
        chrom_len: u64,
        decimals: usize,
    }

    fn fixed_size_test_data() -> FixedSizeTestData {
        let bin_size = 50;
        let chrom_len = 500;
        let tile = make_tile("chrSim", 0, 0, 200, 0, 260);
        let peaks = vec![
            make_peak("chrSim", 10, 2.5),
            make_peak("chrSim", 35, 4.0),
            make_peak("chrSim", 70, 6.5),
            make_peak("chrSim", 115, 5.0),
            make_peak("chrSim", 160, 7.5),
        ];
        let windows = build_fixed_windows(
            bin_size,
            chrom_len,
            tile.core_start() as u64,
            tile.core_end() as u64,
        );

        FixedSizeTestData {
            tile,
            peaks,
            windows,
            decimals: 2,
        }
    }

    fn two_tile_test_data() -> MultiTileTestData {
        let bin_size = 60;
        let chrom_len = 360;
        let tiles = vec![
            make_tile("chrSim", 0, 0, 180, 0, 210),
            make_tile("chrSim", 1, 180, 360, 150, 390),
        ];

        let peaks_by_tile = vec![
            vec![
                make_peak("chrSim", 15, 2.0),
                make_peak("chrSim", 65, 4.5),
                make_peak("chrSim", 145, 6.0),
                make_peak("chrSim", 175, 5.5),
            ],
            vec![
                make_peak("chrSim", 185, 3.5),
                make_peak("chrSim", 225, 4.0),
                make_peak("chrSim", 245, 7.0),
                make_peak("chrSim", 305, 8.0),
            ],
        ];

        let all_windows = build_fixed_windows(bin_size, chrom_len, 0, chrom_len);

        MultiTileTestData {
            tiles,
            peaks_by_tile,
            all_windows,
            bin_size,
            chrom_len,
            decimals: 2,
        }
    }

    fn buffered_unique_rows(
        tile: &Tile,
        windows: &[IndexedInterval<u64>],
        peak_path: &Path,
        decimals: usize,
    ) -> String {
        let mut accumulator =
            WindowAccumulator::new(PeaksWindowAction::OnlyIncludeThesePositionsUnique, decimals);
        accumulator.reset_for_chromosome(tile.chr.clone());
        let mut next_idx = 0usize;
        accumulator.add_windows_for_tile(
            windows,
            &mut next_idx,
            tile.core_start() as u64,
            tile.core_end() as u64,
        );
        stream_tile_peaks(peak_path, |peak| {
            accumulator.push_peak(&peak);
            Ok(())
        })
        .expect("stream peaks for buffered path");
        let mut out = Vec::new();
        accumulator
            .flush_all(&mut out)
            .expect("flush buffered unique windows");
        String::from_utf8(out).expect("buffered unique rows valid utf8")
    }

    fn aligned_unique_rows(chr: &str, peak_path: &Path, decimals: usize) -> String {
        let best = WindowOutputWriter::collect_aligned_unique_peaks(peak_path)
            .expect("collect aligned unique peaks");
        let mut out = String::new();
        for (pos, height) in best {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                chr,
                pos,
                pos + 1,
                format_number(height, decimals)
            ));
        }
        out
    }

    fn buffered_unique_rows_multi(
        tiles: &[Tile],
        windows: &[IndexedInterval<u64>],
        peak_paths: &[PathBuf],
        decimals: usize,
    ) -> String {
        assert_eq!(tiles.len(), peak_paths.len());
        if tiles.is_empty() {
            return String::new();
        }
        let mut accumulator =
            WindowAccumulator::new(PeaksWindowAction::OnlyIncludeThesePositionsUnique, decimals);
        accumulator.reset_for_chromosome(tiles[0].chr.clone());
        let mut next_idx = 0usize;
        let mut out = Vec::new();

        for (tile, path) in tiles.iter().zip(peak_paths.iter()) {
            accumulator.add_windows_for_tile(
                windows,
                &mut next_idx,
                tile.core_start() as u64,
                tile.core_end() as u64,
            );
            stream_tile_peaks(path, |peak| {
                accumulator.push_peak(&peak);
                Ok(())
            })
            .expect("stream peaks for buffered multi unique");
            accumulator
                .flush_completed_windows(tile.core_end() as u64, &mut out)
                .expect("flush completed multi unique windows");
        }
        accumulator
            .flush_all(&mut out)
            .expect("flush remaining multi unique windows");
        String::from_utf8(out).expect("buffered multi unique utf8")
    }

    fn aligned_unique_rows_multi(
        tiles: &[Tile],
        peak_paths: &[PathBuf],
        decimals: usize,
    ) -> String {
        assert_eq!(tiles.len(), peak_paths.len());
        let mut out = String::new();
        for (tile, path) in tiles.iter().zip(peak_paths.iter()) {
            let best = WindowOutputWriter::collect_aligned_unique_peaks(path)
                .expect("collect aligned unique peaks multi");
            for (pos, height) in best {
                out.push_str(&format!(
                    "{}\t{}\t{}\t{}\n",
                    tile.chr,
                    pos,
                    pos + 1,
                    format_number(height, decimals)
                ));
            }
        }
        out
    }

    fn buffered_stats_rows(
        tile: &Tile,
        windows: &[IndexedInterval<u64>],
        peak_path: &Path,
        decimals: usize,
    ) -> String {
        let mut accumulator = WindowAccumulator::new(PeaksWindowAction::Stats, decimals);
        accumulator.reset_for_chromosome(tile.chr.clone());
        let mut next_idx = 0usize;
        accumulator.add_windows_for_tile(
            windows,
            &mut next_idx,
            tile.core_start() as u64,
            tile.core_end() as u64,
        );
        stream_tile_peaks(peak_path, |peak| {
            accumulator.push_peak(&peak);
            Ok(())
        })
        .expect("stream peaks for stats");
        let mut out = Vec::new();
        accumulator
            .flush_completed_windows(tile.core_end() as u64, &mut out)
            .expect("flush completed stat windows");
        accumulator
            .flush_all(&mut out)
            .expect("flush remaining stat windows");
        String::from_utf8(out).expect("buffered stats rows valid utf8")
    }

    fn aligned_stats_rows(
        chr: &str,
        windows: &[IndexedInterval<u64>],
        contributions: &[WindowStatsContribution],
        decimals: usize,
    ) -> String {
        let mut lookup: BTreeMap<u64, &WindowStatsContribution> = BTreeMap::new();
        for contribution in contributions {
            lookup.insert(contribution.window_idx, contribution);
        }
        let mut out = String::new();
        for window in windows {
            let start = window.start();
            let end = window.end();
            let idx = window.idx();
            if let Some(contribution) = lookup.get(&idx) {
                let (avg, median) = stats_distance_summary(
                    contribution.distance_sum,
                    &contribution.distance_histogram,
                );
                out.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    chr,
                    start,
                    end,
                    idx,
                    contribution.count,
                    format_number(avg, decimals),
                    format_number(median, decimals)
                ));
            } else {
                out.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    chr,
                    start,
                    end,
                    idx,
                    0,
                    format_number(f32::NAN, decimals),
                    format_number(f32::NAN, decimals)
                ));
            }
        }
        out
    }

    fn buffered_stats_rows_multi(
        tiles: &[Tile],
        windows: &[IndexedInterval<u64>],
        peak_paths: &[PathBuf],
        decimals: usize,
    ) -> String {
        assert_eq!(tiles.len(), peak_paths.len());
        if tiles.is_empty() {
            return String::new();
        }
        let mut accumulator = WindowAccumulator::new(PeaksWindowAction::Stats, decimals);
        accumulator.reset_for_chromosome(tiles[0].chr.clone());
        let mut next_idx = 0usize;
        let mut out = Vec::new();
        for (tile, path) in tiles.iter().zip(peak_paths.iter()) {
            accumulator.add_windows_for_tile(
                windows,
                &mut next_idx,
                tile.core_start() as u64,
                tile.core_end() as u64,
            );
            stream_tile_peaks(path, |peak| {
                accumulator.push_peak(&peak);
                Ok(())
            })
            .expect("stream peaks for buffered multi stats");
            accumulator
                .flush_completed_windows(tile.core_end() as u64, &mut out)
                .expect("flush completed multi stats windows");
        }
        accumulator
            .flush_all(&mut out)
            .expect("flush remaining multi stats windows");
        String::from_utf8(out).expect("buffered multi stats utf8")
    }

    fn aligned_stats_rows_multi(
        tiles: &[Tile],
        bin_size: u64,
        chrom_len: u64,
        peaks_by_tile: &[Vec<PeakCall>],
        decimals: usize,
    ) -> String {
        assert_eq!(tiles.len(), peaks_by_tile.len());
        let mut out = String::new();
        for (tile, peaks) in tiles.iter().zip(peaks_by_tile.iter()) {
            let windows = build_fixed_windows(
                bin_size,
                chrom_len,
                tile.core_start() as u64,
                tile.core_end() as u64,
            );
            let contributions = compute_window_stats_contributions(windows.as_slice(), peaks);
            out.push_str(&aligned_stats_rows(
                tile.chr.as_str(),
                windows.as_slice(),
                contributions.as_slice(),
                decimals,
            ));
        }
        out
    }

    fn build_fixed_windows(
        bin_size: u64,
        chrom_len: u64,
        tile_start: u64,
        tile_end: u64,
    ) -> Vec<IndexedInterval<u64>> {
        if bin_size == 0 || tile_start >= chrom_len {
            return Vec::new();
        }
        let mut start = (tile_start / bin_size) * bin_size;
        let mut windows = Vec::new();
        while start < tile_end && start < chrom_len {
            let window_start = start;
            let end = (start + bin_size).min(chrom_len);
            let idx = window_start / bin_size;
            windows.push(
                IndexedInterval::new(window_start, end, idx)
                    .expect("test fixed-size windows should be valid"),
            );
            start = start.saturating_add(bin_size);
        }
        windows
    }

    fn write_peaks_file(peaks: &[PeakCall]) -> NamedTempFile {
        let mut temp = NamedTempFile::new().expect("create temp peaks file");
        for peak in peaks {
            writeln!(
                temp,
                "{}\t{}\t{}\t{}\t{}",
                peak.chromosome,
                peak.start(),
                peak.end(),
                peak.peak_position,
                peak.height
            )
            .expect("write peak line");
        }
        temp.flush().expect("flush peak file");
        temp
    }

    fn write_peak_files(peaks_by_tile: &[Vec<PeakCall>]) -> Vec<NamedTempFile> {
        peaks_by_tile
            .iter()
            .map(|peaks| write_peaks_file(peaks))
            .collect()
    }

    fn format_number(value: f32, decimals: usize) -> String {
        if value.is_nan() {
            "NaN".to_string()
        } else {
            format!("{:.*}", decimals, value)
        }
    }
}

#[cfg(test)]
mod tests_peak_signal_processing {
    // MOVE-MODULE-LOCAL: direct private WPS peak signal-processing helper tests.

    use crate::commands::wps_peaks::call_peaks::PeakCall;
    use crate::commands::wps_peaks::wps_peaks::{
        PeakSignalProcessingOptions, compute_window_stats_contributions, peaks_from_wps_values,
    };
    use crate::shared::interval::IndexedInterval;

    fn indexed_windows(entries: &[(u64, u64, u64)]) -> Vec<IndexedInterval<u64>> {
        IndexedInterval::from_tuples(entries).expect("test windows should be valid")
    }

    fn assert_peak(peak: &PeakCall, start: u64, end: u64, height: f32) {
        assert_eq!(peak.start(), start);
        assert_eq!(peak.end(), end);
        assert!(
            (peak.height - height).abs() < 1e-6,
            "expected height {height} got {}",
            peak.height
        );
    }

    #[test]
    fn peaks_from_signal_detects_single_long_run() -> anyhow::Result<()> {
        // Residual WPS has a 55bp plateau starting at index 10, exceeding Snyder's 50bp minimum,
        // so the run should be kept as a single peak once the helper converts residuals into peaks.
        let mut residual = vec![0.0f32; 80];
        for value in residual[10..65].iter_mut() {
            *value = 3.0;
        }
        let opts = PeakSignalProcessingOptions {
            smoothing: false,
            normalization_bp: None,
            min_unmasked: 1,
            min_peak_height: 1.0,
            initial_segment_marker: 0,
        };
        let peaks = peaks_from_wps_values("chrX", 1_000, &residual, None, &opts)?;
        assert_eq!(peaks.len(), 1);
        let peak = &peaks[0];
        assert_peak(peak, 1_010, 1_065, 3.0);
        Ok(())
    }

    #[test]
    fn peaks_from_signal_breaks_runs_on_masked_segments() -> anyhow::Result<()> {
        // Same plateau shape, but we mask a 10bp band (indices 90-99). Snyder requires >=50bp runs,
        // so we extend the positive segments to 10..89 and 100..169 (80bp and 70bp respectively)
        // to stay above the cutoff after the mask splits the trace. Each unmasked run therefore
        // forms its own peak.
        let mut residual = vec![0.0f32; 200];
        for value in residual[10..170].iter_mut() {
            *value = 2.5;
        }
        let mut mask = vec![0u8; residual.len()];
        mask[90..100].fill(1);
        let opts = PeakSignalProcessingOptions {
            smoothing: false,
            normalization_bp: None,
            min_unmasked: 1,
            min_peak_height: 1.0,
            initial_segment_marker: 0,
        };
        let peaks = peaks_from_wps_values("chrY", 500, &residual, Some(&mask), &opts)?;
        assert_eq!(peaks.len(), 2);
        assert_peak(&peaks[0], 510, 590, 2.5);
        assert_peak(&peaks[1], 600, 670, 2.5);
        Ok(())
    }

    #[test]
    fn peaks_from_signal_supports_normalization() -> anyhow::Result<()> {
        // Raw WPS has a 100bp plateau at +5 surrounded by zeros. A 200bp rolling median stays at 0,
        // so residuals remain >0 and the helper should recover one peak covering the plateau.
        let mut wps = vec![0.0f32; 400];
        for value in wps[120..220].iter_mut() {
            *value = 5.0;
        }
        let opts = PeakSignalProcessingOptions {
            smoothing: false,
            normalization_bp: Some(200),
            min_unmasked: 1,
            min_peak_height: 1.0,
            initial_segment_marker: 0,
        };
        let peaks = peaks_from_wps_values("chrZ", 0, &wps, None, &opts)?;
        assert_eq!(peaks.len(), 1);
        let peak = &peaks[0];
        assert_eq!(peak.start(), 120);
        assert_eq!(peak.end(), 220);
        assert!(peak.height > 2.0 && peak.height <= 5.0);
        Ok(())
    }

    #[test]
    fn stats_ignore_distances_across_masked_regions() -> anyhow::Result<()> {
        // Two positive plateaus separated by a masked band emulate two segments inside one window.
        // The stats helper must not report the cross-gap distance because the segment markers differ.
        let mut residual = vec![0.0f32; 200];
        for value in residual[10..70].iter_mut() {
            *value = 2.0;
        }
        for value in residual[130..180].iter_mut() {
            *value = 2.0;
        }
        let mut mask = vec![0u8; residual.len()];
        mask[85..115].fill(1);
        let opts = PeakSignalProcessingOptions {
            smoothing: false,
            normalization_bp: None,
            min_unmasked: 1,
            min_peak_height: 1.0,
            initial_segment_marker: 0,
        };
        let peaks = peaks_from_wps_values("chr1", 0, &residual, Some(&mask), &opts)?;
        assert_eq!(peaks.len(), 2);
        let windows = indexed_windows(&[(0, 200, 0)]);
        let contributions = compute_window_stats_contributions(&windows, &peaks);
        let stats = contributions.first().expect("stats contribution missing");
        assert_eq!(stats.count, 2);
        assert!(stats.distance_histogram.is_empty());
        assert_eq!(stats.distance_sum, 0.0);
        Ok(())
    }

    fn segmented_peak(position: u64, segment: u64) -> PeakCall {
        PeakCall::new(
            "chr1".to_string(),
            position,
            position + 1,
            position,
            1.0,
            segment,
        )
        .expect("test peak should be valid")
    }

    #[test]
    fn peak_call_requires_peak_position_inside_half_open_interval() {
        // Half-open interval semantics are [start, end):
        // - a peak at 99 is before [100,110) and must fail
        // - a peak at 110 lands exactly on the exclusive end and must also fail
        let start_error = PeakCall::new("chr1".to_string(), 100, 110, 99, 1.0, 0)
            .expect_err("peak before interval should fail");
        assert_eq!(
            start_error.to_string(),
            "Peak position 99 must lie inside interval [100, 110)"
        );

        let end_error = PeakCall::new("chr1".to_string(), 100, 110, 110, 1.0, 0)
            .expect_err("peak at exclusive end should fail");
        assert_eq!(
            end_error.to_string(),
            "Peak position 110 must lie inside interval [100, 110)"
        );
    }

    #[test]
    fn stats_remove_single_cross_tile_distance_when_blacklist_hits_boundary() {
        // Tile boundary sits at 2,000bp. First scenario: blacklist begins far upstream so the halo
        // never crosses it, meaning both tiles reuse the same segment marker and the histogram keeps
        // the cross-tile distance 400bp (between the last peak of tile A and first peak of tile B).
        let windows = indexed_windows(&[(0, 4_000, 0)]);
        let peaks_same = vec![
            segmented_peak(1_200, 0),
            segmented_peak(1_500, 0),
            segmented_peak(1_900, 0),
        ];
        let stats_same = compute_window_stats_contributions(&windows, &peaks_same)
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(
            stats_same.distance_histogram.values().copied().sum::<u32>(),
            2,
            "expected both intra- and cross-tile distances when the boundary mask is absent"
        );
        assert_eq!(stats_same.distance_histogram.get(&300), Some(&1));
        assert_eq!(stats_same.distance_histogram.get(&400), Some(&1));

        // Second scenario: move the blacklist right up to the tile edge so tile B sees it in its
        // halo and seeds a new segment marker. Only the intra-tile 300bp distance should remain.
        let mut peaks_masked = peaks_same;
        peaks_masked[2].segment_id = 1_600; // same magnitude as the simulated mask
        let stats_masked = compute_window_stats_contributions(&windows, &peaks_masked)
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(
            stats_masked
                .distance_histogram
                .values()
                .copied()
                .sum::<u32>(),
            1,
            "mask at the boundary should remove exactly one cross-tile distance"
        );
        assert_eq!(stats_masked.distance_histogram.get(&300), Some(&1));
        assert!(
            stats_masked.distance_histogram.get(&400).is_none(),
            "cross-tile 400bp gap must disappear once segments diverge"
        );
    }
}
