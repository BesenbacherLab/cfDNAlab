#![cfg(feature = "cmd_midpoints")]

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use cfdnalab::commands::midpoints::counting_by_group::*;

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
        ProfileGroupsCounts::new(5, 3, vec![20, 50, 100])
    }

    #[test]
    fn new_and_shape() {
        let c = make_counts();
        assert_eq!(c.n_positions(), 5);
        assert_eq!(c.n_groups(), 3);
        assert_eq!(c.n_lengths(), 2);
        assert_eq!(c.counts.len(), 3 * 5 * 2);
        assert_eq!(c.min_fragment_length(), 20);
        assert_eq!(c.max_fragment_length(), 99); // inclusive
    }

    #[test]
    fn index_of_valid_and_bounds() -> Result<()> {
        let c = make_counts();
        // group=1, pos=3, len=20 -> bin 0
        let idx = c.index_of(3, 1, 20)?;
        // idx = group*(P*L) + pos*L + len_bin
        // P=5, L=2 -> 1*(5*2) + 3*2 + 0 = 10 + 6 + 0 = 16
        assert_eq!(idx, 16);

        // Top edge inclusive: 99 ok (bin 1), 100 is out (exclusive)
        assert!(c.index_of(0, 0, 99).is_ok());
        assert!(c.index_of(0, 0, 100).is_err());

        // Below min -> error
        assert!(c.index_of(0, 0, 19).is_err());

        // Position/group out of range -> error
        assert!(c.index_of(5, 0, 20).is_err()); // pos too big
        assert!(c.index_of(0, 3, 20).is_err()); // group too big
        Ok(())
    }

    #[test]
    fn incr_and_get() -> Result<()> {
        let mut c = make_counts();
        // Start zero
        assert_eq!(c.get(0, 0, 20)?, 0.0);
        // Increment integer
        c.incr(0, 0, 20)?;
        assert_eq!(c.get(0, 0, 20)?, 1.0);
        // Increment weighted
        c.incr_weighted(0, 0, 20, 2.5)?;
        approx_eq(c.get(0, 0, 20)?, 3.5, 1e-6);
        // Other bin unaffected
        assert_eq!(c.get(0, 0, 50)?, 0.0);
        Ok(())
    }

    #[test]
    fn incr_errors_on_oob() {
        let mut c = make_counts();
        // Length out of bounds
        assert!(c.incr(0, 0, 19).is_err());
        assert!(c.incr_weighted(0, 0, 100, 1.0).is_err());
        // Position/group out of bounds
        assert!(c.incr(5, 0, 20).is_err());
        assert!(c.incr(0, 3, 20).is_err());
    }

    #[test]
    fn zeroed_like_copies_shape_and_zeros() -> Result<()> {
        let mut c = make_counts();
        c.incr(1, 2, 55)?; // bin 1
        let z = c.zeroed_like();
        assert_eq!(z.counts.len(), c.counts.len());
        assert_eq!(z.n_groups(), c.n_groups());
        assert_eq!(z.n_positions(), c.n_positions());
        assert_eq!(z.n_lengths(), c.n_lengths());
        assert!(z.counts.iter().all(|&v| v == 0.0));
        Ok(())
    }

    #[test]
    fn merge_from_adds_elementwise() -> Result<()> {
        let mut a = make_counts();
        let mut b = make_counts();
        a.incr(2, 1, 49)?; // bin 0
        a.incr_weighted(2, 1, 60, 0.5)?; // bin 1
        b.incr(2, 1, 49)?; // same cells
        b.incr_weighted(2, 1, 60, 1.25)?;
        a.merge_from(&b)?;
        // Check both positions merged
        approx_eq(a.get(2, 1, 49)?, 2.0, 1e-6);
        approx_eq(a.get(2, 1, 60)?, 1.75, 1e-6);
        Ok(())
    }

    #[test]
    fn merge_incompatible_fails() {
        let a = make_counts();
        // Different window size
        let b = ProfileGroupsCounts::new(6, 3, vec![20, 50, 100]);
        assert!(a.clone().merge_from(&b).is_err());
        // Different groups
        let b = ProfileGroupsCounts::new(5, 2, vec![20, 50, 100]);
        assert!(a.clone().merge_from(&b).is_err());
        // Different bins
        let b = ProfileGroupsCounts::new(5, 3, vec![20, 60, 100]);
        assert!(a.clone().merge_from(&b).is_err());
    }

    #[test]
    fn collapse_sums_many() -> Result<()> {
        let mut a = make_counts();
        let mut b = make_counts();
        let mut c = make_counts();
        a.incr(0, 0, 20)?;
        b.incr(0, 0, 20)?;
        c.incr_weighted(0, 0, 20, 3.0)?;
        let total = ProfileGroupsCounts::collapse([&a, &b, &c])?;
        approx_eq(total.get(0, 0, 20)?, 5.0, 1e-6);
        Ok(())
    }

    #[test]
    fn reshape_to_3d_group_len_pos() -> Result<()> {
        let mut c = make_counts(); // G=3, P=5, L=2

        // Fill a few distinct cells we can verify after reshape.
        // (group=2, pos=4, len in bin 0)
        c.incr_weighted(4, 2, 21, 1.5)?;
        // (group=1, pos=0, len in bin 1)
        c.incr(0, 1, 90)?;
        // (group=0, pos=3, len in bin 1)
        c.incr_weighted(3, 0, 60, 2.25)?;

        let m = c.to_3d_group_len_pos();
        assert_eq!(m.len(), 3); // groups
        assert_eq!(m[0].len(), 2); // length bins
        assert_eq!(m[0][0].len(), 5); // positions

        // Check values landed at (group, len_bin, pos)
        approx_eq(m[2][0][4], 1.5, 1e-6); // group 2, bin 0, pos 4
        approx_eq(m[1][1][0], 1.0, 1e-6); // group 1, bin 1, pos 0
        approx_eq(m[0][1][3], 2.25, 1e-6); // group 0, bin 1, pos 3

        // And some zeros elsewhere
        assert_eq!(m[2][1][4], 0.0);
        assert_eq!(m[1][0][0], 0.0);
        Ok(())
    }

    #[test]
    fn ndarray3_view_matches_allocating_copy_for_all_cells() -> Result<()> {
        let mut c = make_counts(); // G=3, P=5, L=2 => flat layout has 30 cells

        // Fill every flat cell with a unique value so any stride mistake shows up immediately.
        // Internal layout is (group, position, length_bin), so the flat index is:
        //   idx = group * (P * L) + position * L + len_bin
        for (idx, value) in c.counts.iter_mut().enumerate() {
            *value = idx as f32 + 0.25;
        }

        let copied = c.to_3d_group_len_pos();
        let viewed = c.view_ndarray3_group_len_pos();

        assert_eq!(viewed.shape(), &[3, 2, 5]);

        // Spot-check a few hand-derived cells before comparing the whole array:
        // - (g=0, len=0, pos=0) -> idx = 0*(5*2) + 0*2 + 0 = 0
        // - (g=1, len=1, pos=3) -> idx = 1*(5*2) + 3*2 + 1 = 17
        // - (g=2, len=0, pos=4) -> idx = 2*(5*2) + 4*2 + 0 = 28
        approx_eq(viewed[(0, 0, 0)], 0.25, 1e-6);
        approx_eq(viewed[(1, 1, 3)], 17.25, 1e-6);
        approx_eq(viewed[(2, 0, 4)], 28.25, 1e-6);

        for group_idx in 0..3 {
            for len_bin in 0..2 {
                for position in 0..5 {
                    approx_eq(
                        viewed[(group_idx, len_bin, position)],
                        copied[group_idx][len_bin][position],
                        1e-6,
                    );
                }
            }
        }

        Ok(())
    }

    #[test]
    fn display_has_shape_info() {
        let c = make_counts();
        let s = format!("{}", c);
        // Basic shape strings are present
        assert!(s.contains("ProfileGroupsCounts("));
        assert!(s.contains("groups:[0..=2]"));
        assert!(s.contains("pos:[0..=4]"));
        // Inclusive max length in the textual form
        assert!(s.contains("len:[20..50...=99]") || s.contains("len:[20..=99]"));
    }
}
