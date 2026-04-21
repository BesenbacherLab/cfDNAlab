use super::{build_summary_prefixes, coverage_sum_and_counts, finalize_value};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::shared::{coverage::Coverage, interval::Interval};

#[test]
fn build_summary_prefixes_without_mask_omits_unmasked_prefixes() {
    // Arrange
    // Use a tiny finalized coverage track:
    //   index: 0   1    2   3
    //   cov:   0, 1.5, 0, 2
    //
    // Hand-derived prefixes:
    // - sum_of_squares_all:
    //     [0,
    //      0^2,
    //      0^2 + 1.5^2,
    //      0^2 + 1.5^2 + 0^2,
    //      0^2 + 1.5^2 + 0^2 + 2^2]
    //   = [0, 0, 2.25, 2.25, 6.25]
    // - nonzero_all:
    //   = [0, 0, 1, 1, 2]
    let mut coverage = Coverage::new(4);
    coverage.finalize_coverage(false);
    coverage
        .coverage_mut()
        .expect("coverage should be available after finalization")
        .copy_from_slice(&[0.0, 1.5, 0.0, 2.0]);

    // Act
    let prefixes = build_summary_prefixes(&coverage).expect("summary prefixes");

    // Assert
    assert_eq!(prefixes.sum_of_squares_all, vec![0.0, 0.0, 2.25, 2.25, 6.25]);
    assert_eq!(prefixes.nonzero_all, vec![0, 0, 1, 1, 2]);
    assert_eq!(prefixes.sum_of_squares_unmasked, None);
    assert_eq!(prefixes.nonzero_unmasked, None);
}

#[test]
fn build_summary_prefixes_with_mask_tracks_all_and_unmasked_prefixes() {
    // Arrange
    // Finalized coverage track:
    //   index: 0  1  2  3
    //   cov:   1, 0, 2, 3
    //
    // Blacklist [1, 3), so allowed positions are indices 0 and 3.
    //
    // Hand-derived prefixes:
    // - sum_of_squares_all:
    //   [0, 1, 1, 5, 14]
    // - nonzero_all:
    //   [0, 1, 1, 2, 3]
    // - sum_of_squares_unmasked:
    //   allowed values are [1, 3], so
    //   [0, 1, 1, 1, 10]
    // - nonzero_unmasked:
    //   [0, 1, 1, 1, 2]
    let mut coverage = Coverage::new(4);
    coverage.finalize_coverage(false);
    coverage
        .coverage_mut()
        .expect("coverage should be available after finalization")
        .copy_from_slice(&[1.0, 0.0, 2.0, 3.0]);
    coverage
        .set_blacklist_mask(&[Interval::new(1, 3).expect("valid blacklist interval")])
        .expect("blacklist mask");

    // Act
    let prefixes = build_summary_prefixes(&coverage).expect("summary prefixes");

    // Assert
    assert_eq!(prefixes.sum_of_squares_all, vec![0.0, 1.0, 1.0, 5.0, 14.0]);
    assert_eq!(prefixes.nonzero_all, vec![0, 1, 1, 2, 3]);
    assert_eq!(
        prefixes.sum_of_squares_unmasked,
        Some(vec![0.0, 1.0, 1.0, 1.0, 10.0])
    );
    assert_eq!(prefixes.nonzero_unmasked, Some(vec![0, 1, 1, 1, 2]));
}

#[test]
fn coverage_sum_and_counts_uses_allowed_prefixes_when_masked_indexes_exist() {
    // Arrange
    // Per-base coverage is [1, 2, 3, 4]
    // Allowed bases are positions 0 and 2, so masked coverage is [1, 0, 3, 0]
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];
    let psum_allowed = [0.0, 1.0, 1.0, 4.0, 4.0];
    let allowed_count_prefix = [0_u32, 1, 1, 2, 2];
    let blacklist_mask = [0_u8, 1, 0, 1];

    // Act
    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
        0,
        4,
        true,
        &psum_all,
        Some(&psum_allowed),
        Some(&allowed_count_prefix),
        Some(&blacklist_mask),
    );

    // Assert
    assert_eq!(sum, 4.0);
    assert_eq!(allowed, 2);
    assert_eq!(blacklisted, 2);
}

#[test]
fn coverage_sum_and_counts_scans_the_mask_when_allowed_count_prefix_is_missing() {
    // Arrange
    // Per-base coverage is still [1, 2, 3, 4]
    // Query [1, 4) -> values [2, 3, 4]
    // Only the middle base is allowed, so:
    //   allowed coverage sum = 3
    //   allowed count = 1
    //   blacklisted count = 2
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];
    let psum_allowed = [0.0, 1.0, 1.0, 4.0, 4.0];
    let blacklist_mask = [0_u8, 1, 0, 1];

    // Act
    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
        1,
        4,
        true,
        &psum_all,
        Some(&psum_allowed),
        None,
        Some(&blacklist_mask),
    );

    // Assert
    assert_eq!(sum, 3.0);
    assert_eq!(allowed, 1);
    assert_eq!(blacklisted, 2);
}

#[test]
fn coverage_sum_and_counts_falls_back_to_full_span_when_mask_support_is_missing() {
    // Arrange
    // Without unmasked prefixes or a blacklist mask, masked mode cannot subtract anything.
    // The helper therefore falls back to treating the whole span as allowed.
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];

    // Act
    let (sum, allowed, blacklisted) =
        coverage_sum_and_counts(1, 3, true, &psum_all, None, None, None);

    // Assert
    assert_eq!(sum, 5.0);
    assert_eq!(allowed, 2);
    assert_eq!(blacklisted, 0);
}

#[test]
fn coverage_sum_and_counts_uses_full_sum_and_span_when_unmasked() {
    // Arrange
    // Per-base coverage is [1, 2, 3, 4]
    // Query [1, 4) -> values [2, 3, 4]
    let psum_all = [0.0, 1.0, 3.0, 6.0, 10.0];
    let psum_allowed = [0.0, 1.0, 1.0, 4.0, 4.0];
    let allowed_count_prefix = [0_u32, 1, 1, 2, 2];
    let blacklist_mask = [0_u8, 1, 0, 1];

    // Act
    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
        1,
        4,
        false,
        &psum_all,
        Some(&psum_allowed),
        Some(&allowed_count_prefix),
        Some(&blacklist_mask),
    );

    // Assert
    assert_eq!(sum, 9.0);
    assert_eq!(allowed, 3);
    assert_eq!(blacklisted, 0);
}

#[test]
fn finalize_value_returns_zero_for_masked_average_with_no_allowed_positions() {
    // Arrange / Act / Assert
    let value = finalize_value(7.5, 0, 100, true, &CoverageWindowAction::Average);
    assert_eq!(value, 0.0);
}

#[test]
fn finalize_value_returns_zero_for_unmasked_average_with_zero_span() {
    // Arrange / Act / Assert
    let value = finalize_value(7.5, 5, 0, false, &CoverageWindowAction::Average);
    assert_eq!(value, 0.0);
}

#[test]
fn finalize_value_returns_sum_for_total_modes_even_when_denominators_are_zero() {
    // Arrange / Act
    let total = finalize_value(7.5, 0, 0, false, &CoverageWindowAction::Total);
    let grouped_total =
        finalize_value(7.5, 0, 0, true, &CoverageWindowAction::TotalOnUniqueBases);

    // Assert
    assert_eq!(total, 7.5);
    assert_eq!(grouped_total, 7.5);
}
