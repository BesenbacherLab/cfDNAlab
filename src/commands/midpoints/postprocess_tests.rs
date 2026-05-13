use super::{
    ProfileLayout, bin_profile, largest_odd_window_that_fits, postprocess_profile,
    smooth_trimmed_profile,
};
use crate::commands::midpoints::smoothing::MidpointSmoothing;
use ndarray::{Array3, array};

#[test]
fn resolve_savgol_fails_with_suggested_window_when_default_does_not_fit() {
    let error = ProfileLayout::resolve(100, 1, MidpointSmoothing::default())
        .expect_err("default 165 bp smoothing should not fit a 100 bp interval");
    let message = error.to_string();

    assert!(
        message.contains("--smoothing savgol=99"),
        "error should suggest the largest odd fitting window, got: {message}"
    );
}

#[test]
fn resolve_explicit_savgol_fails_with_suggested_window_when_it_does_not_fit() {
    let error = ProfileLayout::resolve(100, 1, MidpointSmoothing::SavGol { window_bp: 101 })
        .expect_err("101 bp smoothing should not fit a 100 bp interval");
    let message = error.to_string();

    assert!(
        message.contains("--smoothing savgol=99"),
        "error should suggest the largest odd fitting window, got: {message}"
    );
}

#[test]
fn resolve_savgol_rejects_intervals_shorter_than_7bp() {
    let error = ProfileLayout::resolve(6, 1, MidpointSmoothing::SavGol { window_bp: 5 })
        .expect_err("smoothing should reject output intervals shorter than 7 bp");
    let message = error.to_string();

    assert!(
        message.contains("at least 7 bp"),
        "unexpected short-interval error: {message}"
    );
}

#[test]
fn resolve_savgol_derives_support_radius_as_flank() {
    let layout =
        ProfileLayout::resolve(200, 1, MidpointSmoothing::SavGol { window_bp: 165 })
            .expect("165 bp smoothing should fit a 200 bp interval");

    assert_eq!(layout.smoothing_flank, 82);
    assert_eq!(layout.flanked_length, 364);
}

#[test]
fn largest_odd_window_that_fits_returns_none_below_minimum_interval_size() {
    assert_eq!(largest_odd_window_that_fits(0), None);
    assert_eq!(largest_odd_window_that_fits(6), None);
    assert_eq!(largest_odd_window_that_fits(7), Some(7));
    assert_eq!(largest_odd_window_that_fits(8), Some(7));
    assert_eq!(largest_odd_window_that_fits(9), Some(9));
}

#[test]
fn postprocess_identity_returns_no_owned_copy() {
    let profile = array![[[1.0_f32, 2.0, 3.0]]];
    let layout =
        ProfileLayout::resolve(3, 1, MidpointSmoothing::None)
            .expect("unsmoothed layout should resolve");

    let transformed =
        postprocess_profile(profile.view(), layout).expect("identity postprocess should succeed");

    assert!(transformed.is_none());
}

#[test]
fn postprocess_bins_final_partial_bin_by_actual_width() {
    let profile = array![[[1.0_f32, 2.0, 3.0, 4.0, 5.0]]];
    let layout =
        ProfileLayout::resolve(5, 3, MidpointSmoothing::None)
            .expect("binned layout should resolve");

    let transformed = postprocess_profile(profile.view(), layout)
        .expect("binning should succeed")
        .expect("binning should produce an owned output");

    assert_eq!(transformed.shape(), &[1, 1, 2]);
    assert_eq!(transformed[[0, 0, 0]], 2.0);
    assert_eq!(transformed[[0, 0, 1]], 4.5);
}

#[test]
fn bin_profile_averages_each_group_and_length_bin_with_partial_final_bins() {
    // Directly test the binning kernel, not just the postprocess branch that calls it.
    // With bin_size=2, five input positions become three output bins:
    // [0,2), [2,4), and a one-position final bin [4,5).
    let profile = array![
        [[1.0_f32, 3.0, 5.0, 7.0, 11.0], [2.0, 2.0, 6.0, 10.0, 18.0]],
        [[0.0_f32, 4.0, 8.0, 12.0, 20.0], [5.0, 7.0, 9.0, 13.0, 17.0]]
    ];
    let expected = array![
        [[2.0_f32, 6.0, 11.0], [2.0, 8.0, 18.0]],
        [[2.0_f32, 10.0, 20.0], [6.0, 11.0, 17.0]]
    ];

    let binned = bin_profile(profile.view(), 2);

    assert_eq!(binned.shape(), &[2, 2, 3]);
    assert_eq!(binned, expected);
}

#[test]
fn bin_profile_with_bin_size_one_keeps_position_values() {
    // The caller usually avoids allocating on the identity path, but the helper itself should
    // still have simple, predictable behavior for bin_size=1.
    let profile = array![[[1.0_f32, 4.0, 9.0]]];

    let binned = bin_profile(profile.view(), 1);

    assert_eq!(binned.shape(), &[1, 1, 3]);
    assert_eq!(binned, profile);
}

#[test]
fn postprocess_rejects_unsmoothed_layout_with_extra_counted_positions() {
    // Unsmoothed profiles do not have smoothing flanks. If such a layout claims extra counted positions,
    // silently slicing them away would hide a broken caller-side profile layout.
    let profile = array![[[1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0]]];
    let layout = ProfileLayout {
        output_len: 5,
        flanked_length: 6,
        bin_size: 2,
        smoothing_window: None,
        smoothing_flank: 0,
        output_positions: 3,
    };

    let error = postprocess_profile(profile.view(), layout)
        .expect_err("unsmoothed profiles with extra counted positions should fail");
    let message = error.to_string();

    assert!(
        message.contains("unsmoothed midpoint profile expected 5 positions, got 6"),
        "unexpected unsmoothed profile shape error: {message}"
    );
}

#[test]
fn smooth_trimmed_profile_uses_bases_on_both_sides_of_the_center() {
    // For window=7, the smoothing flank is 3 bases.
    // With output_len=1, the only public output base is counted index 3:
    //
    // Counted index:       0   1   2   3   4   5   6
    // Relative position:  -3  -2  -1   0  +1  +2  +3
    // Coefficient * 21:   -2   3   6   7   6   3  -2
    //
    // Put 21 counts one base left of center and 21 counts one base right of center.
    // Each side contributes 21 * (6 / 21) = 6, so the output must be 12.
    let profile = array![[[0.0_f32, 0.0, 21.0, 0.0, 21.0, 0.0, 0.0]]];

    let smoothed =
        smooth_trimmed_profile(profile.view(), 7, 1).expect("centered support should smooth");

    assert_eq!(smoothed.shape(), &[1, 1, 1]);
    assert!((smoothed[[0, 0, 0]] - 12.0).abs() < 1e-6);
}

#[test]
fn smooth_trimmed_profile_centers_each_output_position_after_the_left_flank() {
    // For output_len=3 and window=7, the counted profile has 3 left-flank bases,
    // 3 public output bases, and 3 right-flank bases:
    //
    // Counted index:        0  1  2 | 3  4  5 | 6  7  8
    // Public output index:          | 0  1  2 |
    //
    // A linear ramp makes an off-by-flank error obvious. The smoother should return the value at
    // the center of each fitted window, so the expected output is the counted value at indexes
    // 3, 4, and 5. This test is about coordinate alignment, not about visible smoothing.
    let values: Vec<f32> = (0..9).map(|value| value as f32).collect();
    let profile =
        Array3::from_shape_vec((1, 1, 9), values).expect("profile shape should fit values");

    let transformed =
        smooth_trimmed_profile(profile.view(), 7, 3).expect("linear profile should smooth");

    assert_eq!(transformed.shape(), &[1, 1, 3]);
    assert!((transformed[[0, 0, 0]] - 3.0).abs() < 1e-6);
    assert!((transformed[[0, 0, 1]] - 4.0).abs() < 1e-6);
    assert!((transformed[[0, 0, 2]] - 5.0).abs() < 1e-6);
}

#[test]
fn smooth_trimmed_profile_smooths_narrow_peak_without_removing_local_trend() {
    // This is closer to a real profile shape: a local trend with one narrow high point.
    //
    // Baseline values 10..16 are linear, so the 7 bp order-3 SavGol smoother preserves the
    // baseline center exactly: 13 at counted index 3.
    //
    // Add a 21-count spike at the center. The center coefficient is 7 / 21, so only 7 counts from
    // that spike remain in the smoothed center. Expected output: 13 + 7 = 20.
    //
    // Raw center count is 34, so this checks the practical behavior: a narrow peak is reduced
    // while the local trend remains in the result.
    let profile = array![[[10.0_f32, 11.0, 12.0, 34.0, 14.0, 15.0, 16.0]]];

    let smoothed =
        smooth_trimmed_profile(profile.view(), 7, 1).expect("spike on local trend should smooth");

    assert_eq!(smoothed.shape(), &[1, 1, 1]);
    assert!((smoothed[[0, 0, 0]] - 20.0).abs() < 1e-6);
}

#[test]
fn postprocess_applies_savgol_to_trimmed_output_positions() {
    // For the first retained output base, the 7 bp order-3 coefficients are
    // [-2, 3, 6, 7, 6, 3, -2] / 21. With input values 1..7, the weighted sum is:
    // (-2*1 + 3*2 + 6*3 + 7*4 + 6*5 + 3*6 - 2*7) / 21 = 84 / 21 = 4
    let profile = array![[[
        1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0
    ]]];
    let layout = ProfileLayout::resolve(7, 1, MidpointSmoothing::SavGol { window_bp: 7 })
        .expect("7 bp smoothing should fit a 7 bp interval");

    let transformed = postprocess_profile(profile.view(), layout)
        .expect("smoothing should succeed")
        .expect("smoothing should produce an owned output");

    assert_eq!(transformed.shape(), &[1, 1, 7]);
    assert!((transformed[[0, 0, 0]] - 4.0).abs() < 1e-6);
}

#[test]
fn postprocess_smooths_before_final_binning() {
    // The 7 bp SavGol filter preserves the linear input exactly at retained centers.
    // The retained values are therefore 4, 5, 6, 7, 8, 9, 10 before final binning.
    // With bin size 3, the expected averages are 5, 8, and 10.
    let profile = array![[[
        1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0
    ]]];
    let layout = ProfileLayout::resolve(7, 3, MidpointSmoothing::SavGol { window_bp: 7 })
        .expect("7 bp smoothing and 3 bp binning should resolve");

    let transformed = postprocess_profile(profile.view(), layout)
        .expect("smoothing plus binning should succeed")
        .expect("smoothing plus binning should produce an owned output");

    assert_eq!(transformed.shape(), &[1, 1, 3]);
    assert!((transformed[[0, 0, 0]] - 5.0).abs() < 1e-6);
    assert!((transformed[[0, 0, 1]] - 8.0).abs() < 1e-6);
    assert!((transformed[[0, 0, 2]] - 10.0).abs() < 1e-6);
}

#[test]
fn postprocess_savgol_preserves_quadratic_profile_at_retained_centers() {
    // The input positions are x=0..12 with values x^2. A 7 bp SavGol window has radius 3, so the
    // retained output positions are centered at x=3..9. Because order-3 SavGol preserves every
    // quadratic at the center, selected retained values should be 3^2, 5^2, and 9^2.
    let values: Vec<f32> = (0..13).map(|position| (position * position) as f32).collect();
    let profile =
        Array3::from_shape_vec((1, 1, 13), values).expect("profile shape should fit values");
    let layout = ProfileLayout::resolve(7, 1, MidpointSmoothing::SavGol { window_bp: 7 })
        .expect("7 bp smoothing should fit a 7 bp retained profile with flanks");

    let transformed = postprocess_profile(profile.view(), layout)
        .expect("quadratic smoothing should succeed")
        .expect("smoothing should produce an owned output");

    assert_eq!(transformed.shape(), &[1, 1, 7]);
    assert!((transformed[[0, 0, 0]] - 9.0).abs() < 1e-5);
    assert!((transformed[[0, 0, 2]] - 25.0).abs() < 1e-5);
    assert!((transformed[[0, 0, 6]] - 81.0).abs() < 1e-5);
}

#[test]
fn postprocess_smooths_each_group_and_length_bin_independently() {
    // A constant profile is preserved by the SavGol coefficients because they sum to 1.
    // Different constants in each group/length-bin cell should therefore stay separate.
    let profile = array![
        [[
            2.0_f32, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0
        ],
        [
            4.0_f32, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0
        ]],
        [[
            6.0_f32, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0
        ],
        [
            8.0_f32, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0
        ]]
    ];
    let layout = ProfileLayout::resolve(7, 1, MidpointSmoothing::SavGol { window_bp: 7 })
        .expect("7 bp smoothing should fit a 7 bp retained profile with flanks");

    let transformed = postprocess_profile(profile.view(), layout)
        .expect("independent smoothing should succeed")
        .expect("smoothing should produce an owned output");

    assert_eq!(transformed.shape(), &[2, 2, 7]);
    assert!((transformed[[0, 0, 0]] - 2.0).abs() < 1e-6);
    assert!((transformed[[0, 1, 3]] - 4.0).abs() < 1e-6);
    assert!((transformed[[1, 0, 4]] - 6.0).abs() < 1e-6);
    assert!((transformed[[1, 1, 6]] - 8.0).abs() < 1e-6);
}
