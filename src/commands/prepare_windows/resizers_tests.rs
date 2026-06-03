mod tests_resizers {
    use crate::commands::prepare_windows::{
        config::{OobPolicy, PrepareConfig},
        resizers::apply_size_transform,
    };
    use crate::shared::interval::Interval;

    fn base_config() -> PrepareConfig {
        let mut cfg = PrepareConfig::default();
        cfg.oob = OobPolicy::Allow;
        cfg
    }

    fn interval(start: u32, end: u32) -> Option<Interval<u32>> {
        Some(Interval::new(start, end).expect("test interval should be valid"))
    }

    #[test]
    fn resize_with_odd_size_centers_window() {
        let mut cfg = base_config();
        cfg.resize = Some(5);
        // Odd input length and odd target size yield a single centered placement
        let transformed = apply_size_transform(10, 21, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (21 - 10) / 2 = 15, size 5 spans 13-18
        assert_eq!(transformed, interval(13, 18));
    }

    #[test]
    fn resize_with_even_size_centers_window_when_parity_matches() {
        let mut cfg = base_config();
        cfg.resize = Some(6);
        // Even input length and even target size yield a single centered placement
        let transformed = apply_size_transform(10, 22, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (22 - 10) / 2 = 16, size 6 spans 13-19
        assert_eq!(transformed, interval(13, 19));
    }

    #[test]
    fn resize_with_even_input_and_odd_target_picks_left_or_right() {
        let mut cfg = base_config();
        cfg.resize = Some(3);
        // Even input length and odd target size have two equally centered placements
        let transformed = apply_size_transform(10, 16, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (16 - 10) / 2 = 13, size 3 can be 11-14 or 12-15
        let left = interval(11, 14);
        let right = interval(12, 15);
        assert!(transformed == left || transformed == right);

        // TODO: Add examples (different seeds) that shows each outcome (regression tests)
    }

    #[test]
    fn resize_with_odd_input_and_even_target_picks_left_or_right() {
        let mut cfg = base_config();
        cfg.resize = Some(4);
        // Odd input length and even target size have two equally centered placements
        let transformed = apply_size_transform(10, 15, Some(100), &cfg).expect("resize");
        // Midpoint is 10 + (15 - 10) / 2 = 12, size 4 can be 10-14 or 11-15
        let left = interval(10, 14);
        let right = interval(11, 15);
        assert!(transformed == left || transformed == right);

        // TODO: Add examples (different seeds) that shows each outcome (regression tests)
    }

    #[test]
    fn flank_with_trim_clamps_to_chrom_bounds() {
        let mut cfg = base_config();
        cfg.flank = Some(vec![5, 5]);
        cfg.oob = OobPolicy::Trim;
        let transformed = apply_size_transform(3, 5, Some(10), &cfg).expect("trim");
        assert_eq!(transformed, interval(0, 10));
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
        assert_eq!(transformed, interval(5, 15));
    }
}
