use super::*;

#[test]
fn deterministic_fragment_midpoint_reuses_same_tie_break() {
    let first_midpoint = midpoint_random_even_for_fragment("chr1", 100, 10);

    for _attempt in 0..32 {
        let repeated_midpoint = midpoint_random_even_for_fragment("chr1", 100, 10);
        assert_eq!(repeated_midpoint, first_midpoint);
    }

    assert!(
        first_midpoint == 104 || first_midpoint == 105,
        "even-length fragment midpoint must be one of the two center bases"
    );
}

#[test]
fn fragment_midpoint_seed_includes_chromosome_start_and_length() {
    let seed = fragment_midpoint_seed("chr1", 100, 10);

    assert_ne!(fragment_midpoint_seed("chr2", 100, 10), seed);
    assert_ne!(fragment_midpoint_seed("chr1", 101, 10), seed);
    assert_ne!(fragment_midpoint_seed("chr1", 100, 12), seed);
}

#[test]
fn deterministic_fragment_midpoint_is_balanced_for_consecutive_starts() {
    let chromosome = "chr1";
    let fragment_length = 100_u32;
    let fragment_count = 10_000_u32;

    let mut left_count = 0_u32;
    let mut right_count = 0_u32;

    for fragment_start in 0..fragment_count {
        let midpoint =
            midpoint_random_even_for_fragment(chromosome, fragment_start, fragment_length);
        let right_center = fragment_start + fragment_length / 2;
        let left_center = right_center - 1;

        if midpoint == left_center {
            left_count += 1;
        } else if midpoint == right_center {
            right_count += 1;
        } else {
            panic!(
                "midpoint {midpoint} is not one of the two centers for fragment start {fragment_start}"
            );
        }
    }

    let minimum_expected_count = fragment_count * 45 / 100;
    let maximum_expected_count = fragment_count * 55 / 100;
    assert!(
        (minimum_expected_count..=maximum_expected_count).contains(&left_count),
        "left midpoint count {left_count} outside 45-55% band for {fragment_count} fragments, right count {right_count}"
    );
}

#[test]
fn odd_length_fragment_midpoint_ignores_seeded_tie_break() {
    let midpoint = midpoint_random_even_for_fragment("chr1", 100, 11);

    assert_eq!(midpoint, 105);
}
