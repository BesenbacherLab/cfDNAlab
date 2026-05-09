#![cfg(feature = "cmd_midpoints")]

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use cfdnalab::commands::midpoints::counting_by_group::*;
    use cfdnalab::shared::length_axis::LengthAxis;
    use std::sync::Arc;

    // Small helper for approximate comparisons where needed
    fn approx_eq(a: f32, b: f32, eps: f32) {
        assert!(
            (a - b).abs() <= eps,
            "expected ~{b}, got {a} (|Δ|={})",
            (a - b).abs()
        );
    }

    fn make_counts() -> ProfileGroupsCounts {
        // window_size=5, groups=3, length bins [20,50), [50,100)
        ProfileGroupsCounts::new(5, 3, length_axis(vec![20, 50, 100]))
    }

    fn length_axis(edges: Vec<u32>) -> Arc<LengthAxis> {
        Arc::new(LengthAxis::new(edges).expect("test length axis should be valid"))
    }

    #[test]
    fn new_and_shape() {
        let counts = make_counts();
        assert_eq!(counts.n_positions(), 5);
        assert_eq!(counts.n_groups(), 3);
        assert_eq!(counts.n_lengths(), 2);
        assert_eq!(counts.counts.len(), 3 * 2 * 5);
        assert_eq!(counts.as_ndarray1().len(), 3 * 2 * 5);
        assert_eq!(counts.min_fragment_length(), 20);
        assert_eq!(counts.max_fragment_length(), 99);
    }

    #[test]
    fn index_of_valid_and_bounds() -> Result<()> {
        let counts = make_counts();

        // Layout is `(group, length_bin, position)`
        // group=1, bin=0, position=3 -> 1*(2*5) + 0*5 + 3 = 13
        assert_eq!(counts.index_of(3, 1, 20)?, 13);
        assert_eq!(counts.index_of(0, 0, 20)?, 0);
        assert_eq!(counts.index_of(4, 0, 20)?, 4);
        assert_eq!(counts.index_of(0, 0, 50)?, 5);
        assert_eq!(counts.index_of(0, 1, 20)?, 10);

        assert!(counts.index_of(0, 0, 99).is_ok());
        assert!(counts.index_of(0, 0, 100).is_err());
        assert!(counts.index_of(0, 0, 19).is_err());
        assert!(counts.index_of(5, 0, 20).is_err());
        assert!(counts.index_of(0, 3, 20).is_err());
        Ok(())
    }

    #[test]
    fn get_reads_count_at_profile_coordinate() -> Result<()> {
        let mut counts = make_counts();

        assert_eq!(counts.get(0, 0, 20)?, 0.0);

        let target_idx = counts.index_of(3, 2, 55)?;
        counts.counts[target_idx] = 2.75;

        approx_eq(counts.get(3, 2, 55)?, 2.75, 1e-6);
        assert_eq!(counts.get(3, 2, 20)?, 0.0);
        Ok(())
    }

    #[test]
    fn get_rejects_out_of_bounds_coordinates() {
        let counts = make_counts();

        assert!(counts.get(0, 0, 19).is_err());
        assert!(counts.get(0, 0, 100).is_err());
        assert!(counts.get(5, 0, 20).is_err());
        assert!(counts.get(0, 3, 20).is_err());
    }

    #[test]
    fn ndarray3_view_exposes_group_length_position_layout() -> Result<()> {
        let mut counts = make_counts();

        let first_idx = counts.index_of(4, 2, 21)?;
        let second_idx = counts.index_of(0, 1, 90)?;
        let third_idx = counts.index_of(3, 0, 60)?;
        counts.counts[first_idx] = 1.5;
        counts.counts[second_idx] = 1.0;
        counts.counts[third_idx] = 2.25;

        let viewed = counts.view_ndarray3_group_len_pos();
        assert_eq!(viewed.shape(), &[3, 2, 5]);

        approx_eq(viewed[(2, 0, 4)], 1.5, 1e-6);
        approx_eq(viewed[(1, 1, 0)], 1.0, 1e-6);
        approx_eq(viewed[(0, 1, 3)], 2.25, 1e-6);

        assert_eq!(viewed[(2, 1, 4)], 0.0);
        assert_eq!(viewed[(1, 0, 0)], 0.0);
        Ok(())
    }

    #[test]
    fn ndarray3_view_matches_flat_index_formula_for_all_coordinates() {
        let mut counts = make_counts();

        // A unique value per flat index makes axis-order mistakes visible
        for (flat_idx, value) in counts.counts.iter_mut().enumerate() {
            *value = flat_idx as f32 + 0.25;
        }

        let viewed = counts.view_ndarray3_group_len_pos();

        assert_eq!(viewed.shape(), &[3, 2, 5]);

        let group_stride: usize = 2 * 5;
        let length_bin_stride: usize = 5;

        // The test fills `counts` by flat index, so the viewed value should be `flat_idx + 0.25`
        let origin_idx = 0 * group_stride + 0 * length_bin_stride + 0;
        let group_1_bin_1_pos_3_idx = group_stride + length_bin_stride + 3;
        let group_2_bin_0_pos_4_idx = 2 * group_stride + 0 * length_bin_stride + 4;

        approx_eq(viewed[(0, 0, 0)], origin_idx as f32 + 0.25, 1e-6);
        approx_eq(
            viewed[(1, 1, 3)],
            group_1_bin_1_pos_3_idx as f32 + 0.25,
            1e-6,
        );
        approx_eq(
            viewed[(2, 0, 4)],
            group_2_bin_0_pos_4_idx as f32 + 0.25,
            1e-6,
        );

        for group_idx in 0..3 {
            for length_bin_idx in 0..2 {
                for position in 0..5 {
                    let flat_idx =
                        group_idx * group_stride + length_bin_idx * length_bin_stride + position;
                    approx_eq(
                        viewed[(group_idx, length_bin_idx, position)],
                        flat_idx as f32 + 0.25,
                        1e-6,
                    );
                }
            }
        }
    }

    #[test]
    fn display_has_shape_info() {
        let counts = make_counts();
        let display_text = format!("{}", counts);

        assert!(display_text.contains("ProfileGroupsCounts("));
        assert!(display_text.contains("groups:[0..=2]"));
        assert!(display_text.contains("pos:[0..=4]"));
        assert!(display_text.contains("len:[20..50...=99]"));
    }
}
