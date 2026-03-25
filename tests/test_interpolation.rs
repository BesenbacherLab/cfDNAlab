use cfdnalab::commands::gc_bias::interpolation::{
    enforce_monotonic_segment, fill_unsupported_bins_with_polynomial,
    fill_zero_bins_with_polynomial,
};

fn assert_slice_close(actual: &[f64], expected: &[f64], tol: f64) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "Slice lengths differ: {} vs {}",
        actual.len(),
        expected.len()
    );
    for (idx, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
        assert!(
            (a - e).abs() <= tol,
            "Mismatch at index {}: actual={} expected={} diff={}",
            idx,
            a,
            e,
            (a - e).abs()
        );
    }
}

fn dense_reference_histogram() -> Vec<f64> {
    vec![
        9.0, 0.0, 0.0, 0.0, 0.0, 19.0, 0.0, 0.0, 0.0, 0.0, 57.0, 0.0, 0.0, 0.0, 0.0, 135.0, 0.0,
        0.0, 0.0, 0.0, 243.0, 0.0, 0.0, 0.0, 0.0, 387.0, 0.0, 0.0, 0.0, 0.0, 510.0, 0.0, 0.0, 0.0,
        0.0, 662.0, 0.0, 0.0, 0.0, 0.0, 728.0, 0.0, 0.0, 0.0, 0.0, 846.0, 0.0, 0.0, 0.0, 0.0,
        1018.0, 0.0, 0.0, 0.0, 0.0, 1196.0, 0.0, 0.0, 0.0, 0.0, 1272.0, 0.0, 0.0, 0.0, 0.0, 1152.0,
        0.0, 0.0, 0.0, 0.0, 946.0, 0.0, 0.0, 0.0, 0.0, 629.0, 0.0, 0.0, 0.0, 0.0, 335.0, 0.0, 0.0,
        0.0, 0.0, 168.0, 0.0, 0.0, 0.0, 0.0, 68.0, 0.0, 0.0, 0.0, 0.0, 15.0, 0.0, 0.0, 0.0, 0.0,
        6.0,
    ]
}

mod zero_interpolator_tests {
    use super::*;

    #[test]
    fn should_leave_histogram_unchanged_without_zero_bins() -> anyhow::Result<()> {
        // Human verification status: unverified
        let mut histogram = vec![2.0, 4.0, 8.0];

        fill_zero_bins_with_polynomial(&mut histogram, 1, 2, 2)?;

        assert_eq!(histogram, vec![2.0, 4.0, 8.0]);

        Ok(())
    }

    #[test]
    fn should_linearly_interpolate_between_two_anchors() -> anyhow::Result<()> {
        // Human verification status: unverified
        // One zero run between two anchor bins; using degree-one interpolation makes
        // manual expectations straightforward: each bin steps by +2.
        let mut histogram = vec![2.0, 0.0, 0.0, 8.0];

        fill_zero_bins_with_polynomial(&mut histogram, 1, 2, 1)?;

        assert_slice_close(&histogram, &[2.0, 4.0, 6.0, 8.0], 1e-9);

        Ok(())
    }

    #[test]
    fn edge_run_with_single_anchor_is_left_unchanged() -> anyhow::Result<()> {
        // Human verification status: unverified
        // Run at the left edge has only one real anchor, so we skip interpolation
        // instead of inventing a slope with fabricated padding.
        let mut histogram = vec![0.0, 0.0, 0.0, 3.0];

        fill_zero_bins_with_polynomial(&mut histogram, 1, 2, 2)?;

        assert_slice_close(&histogram, &[0.0, 0.0, 0.0, 3.0], 1e-9);

        Ok(())
    }

    #[test]
    fn should_not_overshoot_when_zero_run_touches_edge() -> anyhow::Result<()> {
        // Human verification status: unverified
        // Sequence starts with a long zero run followed by a small positive anchor.
        // Clamp logic should prevent us from creating inflated counts before the anchor.
        let mut histogram = vec![0.0, 0.0, 0.0, 0.0, 5.0, 7.0, 11.0, 0.0, 0.0, 0.0, 0.0, 15.0];

        fill_zero_bins_with_polynomial(&mut histogram, 2, 3, 3)?;

        assert!(
            histogram[..4]
                .iter()
                .enumerate()
                .all(|(idx, &value)| value >= 0.0 && value <= 5.0 * (idx as f64 + 1.0)),
            "Edge interpolation should stay within progressively larger bounds"
        );

        Ok(())
    }

    #[test]
    fn should_fill_all_zero_bins_in_reference_example() -> anyhow::Result<()> {
        // Human verification status: unverified
        // Full-length GC count profile from the prompt; interpolation should clear
        // out zero runs while preserving the original anchor bins.
        let mut histogram = dense_reference_histogram();

        fill_zero_bins_with_polynomial(&mut histogram, 2, 3, 3)?;

        // Regression test: store the full interpolated profile to guard against behavioural changes.
        // The expected values were recorded manually from the current implementation.
        let expected = vec![
            9.0,
            9.0,
            9.892320131239746,
            12.630867510350752,
            16.369729318022028,
            19.0,
            19.991333489571108,
            27.054950394500395,
            35.51499856782018,
            45.371478009530456,
            57.0,
            70.52397195555398,
            84.35700466943274,
            99.55817708395297,
            116.12748919911468,
            135.0,
            155.76644866507124,
            175.7483149149669,
            196.80060881988751,
            218.923330379833,
            243.0,
            274.6668036679938,
            300.00727620484565,
            325.88139093699334,
            352.2891478644367,
            387.0,
            411.22663220408657,
            436.82045734908934,
            462.2662996002939,
            487.56415895769993,
            510.0,
            548.1148301474614,
            574.3349600135143,
            599.8823742701591,
            624.7570729173956,
            662.0,
            663.0083903830929,
            683.6670234515229,
            704.266932357445,
            724.8081171008586,
            728.0,
            749.171893366406,
            773.0300656671791,
            798.0205686878347,
            824.1434024283727,
            846.0,
            887.2180598358951,
            916.393499156519,
            946.1848329889024,
            976.592061333045,
            1018.0,
            1071.9311608012918,
            1100.6772721796237,
            1127.2591212628795,
            1151.6767080510617,
            1196.0,
            1227.6873623151669,
            1241.59938924398,
            1250.7090070494978,
            1255.0162157317172,
            1272.0,
            1249.5940159127222,
            1236.3610107692657,
            1217.854399068343,
            1194.0741808099538,
            1152.0,
            1115.4926048028265,
            1079.8948121973444,
            1040.4844607892883,
            997.2615505786544,
            946.0,
            868.6744262499292,
            816.1776740386258,
            762.6059182635336,
            707.9591589246529,
            629.0,
            556.5569655782147,
            503.59266398329055,
            452.5982127091311,
            403.5736117557335,
            335.0,
            300.53201001062007,
            261.18140414522895,
            225.05235209300372,
            192.1448538539462,
            168.0,
            145.52157773323597,
            120.75489274328902,
            98.57088059737907,
            78.96954129549522,
            68.0,
            53.669690847083984,
            40.844285670364116,
            30.020030859139297,
            21.196926413414985,
            15.0,
            13.079555892631106,
            8.58987930490457,
            6.0,
            6.0,
            6.0,
        ];
        assert_slice_close(&histogram, &expected, 1e-9);

        Ok(())
    }
}

mod unsupported_interpolator_tests {
    use super::*;

    #[test]
    fn unsupported_interp_leaves_supported_bins_unchanged() -> anyhow::Result<()> {
        // Human verification status: unverified
        let mut histogram = vec![2.0, 4.0, 8.0];
        let mut mask = vec![true, true, true];

        fill_unsupported_bins_with_polynomial(&mut histogram, &mut mask, 1, 2, 2, true)?;

        assert_eq!(histogram, vec![2.0, 4.0, 8.0]);
        assert_eq!(mask, vec![true, true, true]);

        Ok(())
    }

    #[test]
    fn unsupported_interp_fills_gaps_between_supported_bins() -> anyhow::Result<()> {
        // Human verification status: unverified
        let mut histogram = vec![2.0, 0.0, 0.0, 8.0];
        let mut mask = vec![true, false, false, true];

        fill_unsupported_bins_with_polynomial(&mut histogram, &mut mask, 1, 2, 1, true)?;

        assert_slice_close(&histogram, &[2.0, 4.0, 6.0, 8.0], 1e-9);
        assert!(mask.iter().all(|supported| *supported));

        Ok(())
    }

    #[test]
    fn unsupported_interp_edge_run_with_single_anchor_is_left_unchanged() -> anyhow::Result<()> {
        // Human verification status: unverified
        // Run at the left edge has only one supported anchor, interpolation is
        // skipped to avoid fabricating a slope from padded zeros.
        let mut histogram = vec![0.0, 0.0, 0.0, 3.0];
        let mut mask = vec![false, false, false, true];

        fill_unsupported_bins_with_polynomial(&mut histogram, &mut mask, 1, 2, 2, true)?;

        assert_slice_close(&histogram, &[0.0, 0.0, 0.0, 3.0], 1e-9);
        assert_eq!(mask, vec![false, false, false, true]);

        Ok(())
    }

    #[test]
    fn unsupported_interp_requires_two_real_anchors() -> anyhow::Result<()> {
        // Human verification status: unverified
        // With only one real anchor the solver bails out early, so the unsupported
        // bins remain zeros even though padding is requested.
        let mut histogram = vec![0.0, 0.0, 5.0];
        let mut mask = vec![false, false, true];

        fill_unsupported_bins_with_polynomial(&mut histogram, &mut mask, 1, 2, 1, true)?;

        assert_slice_close(&histogram, &[0.0, 0.0, 5.0], 1e-9);
        assert_eq!(mask, vec![false, false, true]);

        Ok(())
    }

    #[test]
    fn unsupported_interp_interpolates_when_both_sides_supported() -> anyhow::Result<()> {
        // Human verification status: unverified
        // Both ends have genuine anchors (mask entries true), so interpolation
        // fills the interior run with values between 5 and 10.
        let mut histogram = vec![5.0, 0.0, 0.0, 10.0];
        let mut mask = vec![true, false, false, true];

        fill_unsupported_bins_with_polynomial(&mut histogram, &mut mask, 1, 2, 1, true)?;

        assert_slice_close(
            &histogram,
            &[5.0, 6.666666666666667, 8.333333333333334, 10.0],
            1e-9,
        );
        assert_eq!(mask, vec![true, true, true, true]);

        Ok(())
    }

    #[test]
    fn unsupported_interp_skips_when_insufficient_anchors() -> anyhow::Result<()> {
        // Human verification status: unverified
        let mut histogram = vec![0.0, 0.0];
        let mut mask = vec![false, false];

        fill_unsupported_bins_with_polynomial(&mut histogram, &mut mask, 1, 2, 1, true)?;

        assert_eq!(histogram, vec![0.0, 0.0]);
        assert_eq!(mask, vec![false, false]);

        Ok(())
    }

    #[test]
    fn unsupported_interp_matches_zero_interpolator_on_dense_example() -> anyhow::Result<()> {
        // Human verification status: unverified
        let mut zero_histogram = dense_reference_histogram();
        fill_zero_bins_with_polynomial(&mut zero_histogram, 2, 3, 3)?;

        let mut unsupported_histogram = dense_reference_histogram();
        let mut mask: Vec<bool> = unsupported_histogram.iter().map(|v| *v > 0.0).collect();
        fill_unsupported_bins_with_polynomial(
            &mut unsupported_histogram,
            &mut mask,
            2,
            3,
            3,
            true,
        )?;

        assert_slice_close(&unsupported_histogram, &zero_histogram, 1e-9);
        assert!(mask.iter().all(|supported| *supported));

        Ok(())
    }
}

#[test]
fn enforces_non_decreasing_segments() {
    // Human verification status: unverified
    let mut segment = vec![1.0, 0.5, 0.6, 1.2];
    enforce_monotonic_segment(&mut segment, 1.0, 2.0);
    assert_eq!(segment, vec![1.0, 1.0, 1.0, 1.2]);
}

#[test]
fn enforces_non_increasing_segments() {
    // Human verification status: unverified
    let mut segment = vec![5.0, 6.0, 4.0, 3.0];
    enforce_monotonic_segment(&mut segment, 5.0, 1.0);
    assert_eq!(segment, vec![5.0, 5.0, 4.0, 3.0]);
}
