use super::create_support_mask_threshold_per_mb;
use anyhow::Result;
use ndarray::array;

#[test]
fn support_threshold_per_mb_steps_at_hundred_million_positions() -> Result<()> {
    // Arrange:
    // Use exactly 1 Mb of valid ACGT positions so the absolute threshold equals `threshold_per_mb`.
    // The row values are chosen so each threshold step changes exactly one additional support bit.
    let counts = array![[0.5_f64, 1.5_f64, 2.5_f64, 3.5_f64]];
    let num_acgt_positions = 1_000_000_u64;

    let scenarios = [
        (99_999_999_usize, 1_usize, vec![false, true, true, true]),
        (100_000_000_usize, 2_usize, vec![false, false, true, true]),
        (200_000_000_usize, 3_usize, vec![false, false, false, true]),
    ];

    for (n_positions, expected_threshold_per_mb, expected_mask_row) in scenarios {
        // Act:
        // The command computes:
        //   threshold_per_mb = 1 + n_positions / 100_000_000
        // with integer division.
        let threshold_per_mb = 1 + n_positions / 100_000_000;
        let mask = create_support_mask_threshold_per_mb(
            std::slice::from_ref(&counts),
            num_acgt_positions,
            threshold_per_mb as f64,
        )
        .expect("support mask should be created");

        // Assert:
        // Because num_acgt_positions is exactly 1 Mb, the usable-count threshold is:
        //   threshold = 1.0 * threshold_per_mb
        // So the support mask should flip at the exact crossover boundaries.
        assert_eq!(threshold_per_mb, expected_threshold_per_mb);
        assert_eq!(mask.row(0).to_vec(), expected_mask_row);
    }

    Ok(())
}
