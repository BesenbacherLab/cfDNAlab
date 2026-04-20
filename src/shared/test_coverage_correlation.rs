use crate::shared::interval::{TouchingMergePolicy, merge_intervals};

fn build_mask_from_intervals(length: usize, intervals: &[Interval<u32>]) -> Vec<u8> {
    // Convert checked half-open intervals into the explicit positional 0/1 mask.
    //
    // This is the direct binary representation used by the standard Pearson calculation:
    // - 1 means the base is inside the masked group
    // - 0 means the base is outside the masked group
    let mut mask = vec![0u8; length];
    for interval in intervals {
        for pos in interval.start() as usize..interval.end() as usize {
            mask[pos] = 1;
        }
    }
    mask
}

fn build_simulated_coverage(mask: &[u8]) -> Vec<f64> {
    // Build a deterministic positional coverage track with a real signal.
    //
    // Requirements for this fixture:
    // - not random, so the proof stays reproducible
    // - enough positional variation that Pearson is not dominated by a flat baseline
    // - slightly lower coverage inside the mask, so the true correlation is clearly negative
    let mut coverage = Vec::with_capacity(mask.len());

    for pos in 0..mask.len() {
        // Deterministic structured signal over a few hundred positions.
        // Bases inside the mask are shifted downward so the final Pearson R is negative
        // and clearly away from zero.
        let baseline = 18.0
            + (pos % 19) as f64 * 0.45
            + ((pos * 7) % 11) as f64 * 0.18
            + if pos % 31 == 0 { 0.8 } else { 0.0 };
        let masked_shift = if mask[pos] == 1 { 3.4 } else { 0.0 };
        coverage.push(baseline - masked_shift);
    }

    coverage
}

fn pearson_r_from_positional_mask(coverage: &[f64], mask: &[u8]) -> f64 {
    // Ordinary positional Pearson correlation over the explicit vectors.
    //
    // Formula:
    //   R = sum_i[(x_i - mean_x) * (y_i - mean_y)] /
    //       sqrt(sum_i[(x_i - mean_x)^2] * sum_i[(y_i - mean_y)^2])
    //
    // Here:
    // - x_i is the coverage at base i
    // - y_i is the binary 0/1 mask value at base i
    //
    // This is the reference calculation that the interval-derived shortcut must match.
    let n = coverage.len() as f64;

    let mean_x: f64 = coverage.iter().sum::<f64>() / n;
    let mean_y: f64 = mask.iter().map(|&value| value as f64).sum::<f64>() / n;

    let mut covariance_numer = 0.0_f64;
    let mut variance_x_numer = 0.0_f64;
    let mut variance_y_numer = 0.0_f64;

    for (coverage_value, &mask_value) in coverage.iter().zip(mask.iter()) {
        let centered_x = coverage_value - mean_x;
        let centered_y = mask_value as f64 - mean_y;
        covariance_numer += centered_x * centered_y;
        variance_x_numer += centered_x * centered_x;
        variance_y_numer += centered_y * centered_y;
    }

    covariance_numer / (variance_x_numer * variance_y_numer).sqrt()
}

fn pearson_r_from_intervals(coverage: &[f64], merged_intervals: &[Interval<u32>]) -> f64 {
    // Interval-derived Pearson correlation using only global coverage moments plus
    // grouped interval sums after collapsing to unique bases.
    //
    // Let:
    // - n  = total number of positions
    // - S  = sum_i x_i
    // - Q  = sum_i x_i^2
    // - n1 = sum_i y_i, where y_i is the binary 0/1 group mask
    // - S1 = sum_i x_i * y_i
    //
    // Then Pearson correlation between x and the binary mask y can be written as:
    //
    //   R = (n*S1 - S*n1) / sqrt((n*Q - S^2) * (n*n1 - n1^2))
    //
    // This helper computes exactly that form, with `merged_intervals` representing the
    // unique bases of the group.
    let n = coverage.len() as f64;
    let sum_x: f64 = coverage.iter().sum();
    let sum_x2: f64 = coverage.iter().map(|value| value * value).sum();

    let mut n1 = 0.0_f64;
    let mut sum_xy = 0.0_f64;
    for interval in merged_intervals {
        // Because the intervals are already merged, every covered base contributes once to:
        // - n1 = number of bases inside the binary mask
        // - S1 = sum of coverage values inside the binary mask
        n1 += interval.len() as f64;
        for pos in interval.start() as usize..interval.end() as usize {
            sum_xy += coverage[pos];
        }
    }

    let numerator = n * sum_xy - sum_x * n1;
    let denominator = ((n * sum_x2 - sum_x * sum_x) * (n * n1 - n1 * n1)).sqrt();
    numerator / denominator
}

fn site_weighted_r_from_raw_intervals(coverage: &[f64], raw_intervals: &[Interval<u32>]) -> f64 {
    // The same aggregated formula as above, but applied directly to the raw grouped intervals
    // without collapsing overlaps first.
    //
    // This is intentionally a different quantity:
    // - overlap bases contribute multiple times to n1 and S1
    // - the result is therefore site-weighted, not binary-mask based
    //
    // We keep this helper to prove the negative case: raw overlapping grouped intervals should
    // not reproduce the binary-mask Pearson R.
    let n = coverage.len() as f64;
    let sum_x: f64 = coverage.iter().sum();
    let sum_x2: f64 = coverage.iter().map(|value| value * value).sum();

    let mut n1 = 0.0_f64;
    let mut sum_xy = 0.0_f64;
    for interval in raw_intervals {
        n1 += interval.len() as f64;
        for pos in interval.start() as usize..interval.end() as usize {
            sum_xy += coverage[pos];
        }
    }

    let numerator = n * sum_xy - sum_x * n1;
    let denominator = ((n * sum_x2 - sum_x * sum_x) * (n * n1 - n1 * n1)).sqrt();
    numerator / denominator
}

#[test]
fn given_overlapping_mask_intervals_when_merged_then_interval_formula_matches_positional_pearson() {
    // Human verification status: unverified
    //
    // Arrange
    // -------
    // Start from interval masks, not from a prebuilt 0/1 vector.
    // The input intentionally includes both overlapping and touching intervals so the
    // binary-mask interpretation requires merging before counting unique bases:
    //   [24, 72) overlaps [58, 91)
    //   [58, 91) touches  [91, 118)
    //   [210, 244) overlaps [236, 279)
    //
    // After converting those intervals to a 0/1 mask, each base can only be 0 or 1.
    //
    // So this proof compares two mathematically equivalent paths:
    //
    // 1. Positional path
    //    - convert merged intervals to the explicit 0/1 mask vector
    //    - compute ordinary positional Pearson over `(coverage_i, mask_i)`
    //
    // 2. Interval path
    //    - keep only:
    //      n  = total positions
    //      S  = sum coverage globally
    //      Q  = sum squared coverage globally
    //      n1 = unique bases inside the merged intervals
    //      S1 = sum coverage inside the merged intervals
    //    - compute:
    //      R = (n*S1 - S*n1) / sqrt((n*Q - S^2) * (n*n1 - n1^2))
    //
    // If the derivation is correct, the two values must match to floating-point precision.
    let raw_intervals = Interval::from_tuples(&[
        (24_u32, 72_u32),
        (58, 91),
        (91, 118),
        (150, 186),
        (210, 244),
        (236, 279),
        (320, 348),
    ])
    .expect("fixture intervals must be valid");
    let merged_intervals = merge_intervals(raw_intervals.clone(), TouchingMergePolicy::MergeTouching);
    let mask = build_mask_from_intervals(384, &merged_intervals);
    let coverage = build_simulated_coverage(mask.as_slice());

    // Act
    // ---
    let positional_r = pearson_r_from_positional_mask(coverage.as_slice(), mask.as_slice());
    let interval_r = pearson_r_from_intervals(coverage.as_slice(), merged_intervals.as_slice());

    // Assert
    // ------
    // The coverage was shifted down inside the mask, so the signal should be clearly negative.
    // Then the interval-derived formula must recover the exact same R.
    assert!(
        positional_r < -0.40,
        "expected a clear negative signal from the simulated inside-mask coverage dip, got {positional_r}"
    );
    assert!(
        (positional_r - interval_r).abs() < 1.0e-12,
        "merged interval formula should equal the direct positional Pearson R; positional={positional_r}, interval={interval_r}"
    );
}

#[test]
fn given_overlapping_raw_intervals_when_not_merged_then_site_weighted_formula_differs_from_binary_mask_pearson() {
    // Human verification status: unverified
    //
    // Arrange
    // -------
    // Use the same raw interval fixture as above, but compare:
    // 1) direct Pearson against the 0/1 mask built from unique covered bases
    // 2) site-weighted Pearson that counts overlaps in the raw interval list multiple times
    //
    // These are different questions, so they should not give the same result once overlaps exist:
    //
    // - the positional 0/1 mask asks "is this base in the group or not?"
    // - the raw grouped intervals ask "how much site-weighted membership did this base receive?"
    //
    // A mismatch here is the expected result and demonstrates why overlap collapsing is required
    // before interpreting the grouped intervals as a binary mask.
    let raw_intervals = Interval::from_tuples(&[
        (24_u32, 72_u32),
        (58, 91),
        (91, 118),
        (150, 186),
        (210, 244),
        (236, 279),
        (320, 348),
    ])
    .expect("fixture intervals must be valid");
    let merged_intervals = merge_intervals(raw_intervals.clone(), TouchingMergePolicy::MergeTouching);
    let mask = build_mask_from_intervals(384, &merged_intervals);
    let coverage = build_simulated_coverage(mask.as_slice());

    // Act
    // ---
    let positional_r = pearson_r_from_positional_mask(coverage.as_slice(), mask.as_slice());
    let site_weighted_r =
        site_weighted_r_from_raw_intervals(coverage.as_slice(), raw_intervals.as_slice());

    // Assert
    // ------
    // The site-weighted version double-counts overlap bases, so it should not match the binary
    // mask correlation that the merged unique-base intervals represent.
    assert!(
        (positional_r - site_weighted_r).abs() > 1.0e-3,
        "raw overlapping intervals should not reproduce the binary-mask Pearson R without merging; positional={positional_r}, site_weighted={site_weighted_r}"
    );
}
