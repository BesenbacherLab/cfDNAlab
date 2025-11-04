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
        for (original, filtered) in values.iter().zip(smoothed.iter()) {
            assert!(
                approx_eq(*original, *filtered, 1e-4),
                "expected {original} ~= {filtered}"
            );
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
        assert!(approx_eq(smoothed[27], values[27], 1e-4));
        assert!(approx_eq(smoothed[33], values[33], 1e-4));
    }

    #[test]
    fn normalize_subtracts_sliding_median() {
        let numerator = vec![1.0, 2.0, 50.0, 4.0, 5.0, 6.0];
        let baseline = numerator.clone();
        let mask = vec![0u8; numerator.len()];

        let normalized = normalize_wps(&numerator, &baseline, Some(&mask), 5, 1, 3);

        let expected = vec![-1.0, -1.0, 46.0, -1.0, -0.5, 1.0];
        for (idx, (observed, exp)) in normalized.iter().zip(expected.iter()).enumerate() {
            if exp.is_nan() {
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

    const SG_WINDOW_SIZE: u32 = 21;
    const EPS: f32 = 1e-6;

    fn approx_eq(a: f32, b: f32) -> bool {
        (a - b).abs() <= EPS
    }

    #[test]
    fn build_left_edge_window_reflects_prefix() {
        let edge_slice = vec![3.0, 4.0, 5.0];
        let window = build_left_edge_window(&edge_slice);
        assert_eq!(window.len(), SG_WINDOW_SIZE);
        assert!(approx_eq(window[0], 2.0));
        assert!(approx_eq(window[1], 1.0));
        for value in &window[2..SG_WINDOW_SIZE - edge_slice.len()] {
            assert!(approx_eq(*value, 3.0));
        }
        assert_eq!(
            &window[SG_WINDOW_SIZE - edge_slice.len()..],
            edge_slice.as_slice()
        );
    }

    #[test]
    fn build_right_edge_window_reflects_suffix() {
        let edge_slice = vec![1.0, 2.0, 3.0];
        let window = build_right_edge_window(&edge_slice);
        assert_eq!(window.len(), SG_WINDOW_SIZE);
        assert_eq!(&window[..edge_slice.len()], edge_slice.as_slice());
        assert!(approx_eq(window[edge_slice.len()], 4.0));
        assert!(approx_eq(window[edge_slice.len() + 1], 5.0));
        for value in &window[edge_slice.len() + 2..] {
            assert!(approx_eq(*value, 3.0));
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
