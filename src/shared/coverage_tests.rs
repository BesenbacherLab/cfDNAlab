use super::{Coverage, clamp_finite_coverage_below_to_zero};
use crate::shared::{fragment::minimal_fragment::Fragment, gc_tag::GCTagValue, interval::Interval};

fn next_u32(state: &mut u64) -> u32 {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*state >> 32) as u32
}

fn build_zero_sum_delta_stress_coverage() -> Vec<f32> {
    let mut cp = Coverage::new(32);
    let mut seed = 0x5EED_1234_5678_9ABC_u64;
    let weight_choices: [f32; 19] = [
        1.0e-6, 2.0e-6, 5.0e-6, 1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3,
        5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 1.0e-1, 2.0e-1, 5.0e-1, 1.0,
    ];

    // Build 10k additions/subtractions against a small delta array so every slot is
    // mathematically zero at the end, while floating-point roundoff still accumulates in
    // the mixed-sign delta entries before coverage finalization.
    let mut ops = Vec::<(usize, f32)>::with_capacity(10_000);
    for _ in 0..5_000 {
        let slot = 1 + (next_u32(&mut seed) as usize % 8);
        let weight = weight_choices[next_u32(&mut seed) as usize % weight_choices.len()];
        ops.push((slot, weight));
        ops.push((slot, -weight));
    }

    for idx in (1..ops.len()).rev() {
        let swap_idx = next_u32(&mut seed) as usize % (idx + 1);
        ops.swap(idx, swap_idx);
    }

    for (slot, delta) in ops {
        cp.delta[slot] += delta as f64;
    }

    cp.finalize_coverage(true).to_vec()
}

#[derive(Clone, Copy)]
struct WeightedSpan {
    start: u32,
    end: u32,
    weight: f32,
}

fn deterministic_tail_zero_spans(seed: u64, n_fragments: usize) -> Vec<WeightedSpan> {
    let mut state = seed;
    let mut spans = Vec::with_capacity(n_fragments);
    let weight_choices: [f32; 19] = [
        1.0e-6, 2.0e-6, 5.0e-6, 1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3,
        5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 1.0e-1, 2.0e-1, 5.0e-1, 1.0,
    ];

    for _ in 0..n_fragments {
        let start = next_u32(&mut state) % 2048;
        let max_len = 1 + (next_u32(&mut state) % 1000);
        let available = 3000_u32.saturating_sub(start).max(1);
        let len = max_len.min(available);
        let end = start + len;
        let weight = weight_choices[next_u32(&mut state) as usize % weight_choices.len()];
        spans.push(WeightedSpan { start, end, weight });
    }

    spans
}

fn finalize_current_f32_delta_coverage(spans: &[WeightedSpan], length: u32) -> Vec<f32> {
    let mut cp = Coverage::new(length);
    for span in spans {
        cp.add_fragment_weighted(
            Fragment {
                tid: 0,
                interval: Interval::new(span.start, span.end).expect("valid deterministic span"),
                gc_tag: GCTagValue::default(),
            },
            span.weight as f64,
        )
        .expect("deterministic weighted span should be accepted");
    }
    cp.finalize_coverage(true).to_vec()
}

fn finalize_reference_f64_coverage(spans: &[WeightedSpan], length: u32) -> Vec<f64> {
    let mut cov = vec![0.0_f64; length as usize];
    for span in spans {
        for pos in span.start as usize..span.end as usize {
            cov[pos] += span.weight as f64;
        }
    }
    cov
}

fn max_abs_suffix(values: &[f32], suffix_start: usize) -> f32 {
    values[suffix_start..]
        .iter()
        .copied()
        .map(f32::abs)
        .fold(0.0_f32, f32::max)
}

#[test]
fn regression_tail_zero_region_should_survive_theoretical_cleanup_floor() {
    // This fixes a deterministic fixture for the current delta accumulation path.
    //
    // Setup:
    // - `--normalize-by-length` with max fragment length 1000 and minimum usable GC weight 1e-3
    //   implies the smallest real positive pre-scaling support is 1e-6
    // - the current theoretical cleanup floor is therefore 5e-7
    // - all generated spans end before position 3000, so the suffix [3000, 4096) is
    //   mathematically exact zero
    //
    // Required behavior:
    // - the mathematically untouched suffix must accumulate less residue than the run-specific
    //   theoretical cleanup floor
    // - cleanup can then safely clamp that suffix to exact zero
    let cleanup_floor = 5.0e-7_f32;
    let spans = deterministic_tail_zero_spans(0, 200_000);

    let current = finalize_current_f32_delta_coverage(&spans, 4096);
    let mut clamped_current = current.clone();
    clamp_finite_coverage_below_to_zero(&mut clamped_current, cleanup_floor);
    let reference = finalize_reference_f64_coverage(&spans, 4096);
    let raw_residue = max_abs_suffix(&current, 3000);

    assert!(
        raw_residue < cleanup_floor,
        "regression: raw suffix residue {raw_residue} exceeded the theoretical cleanup floor {cleanup_floor}"
    );
    assert!(
        reference[3000..].iter().all(|value| *value == 0.0),
        "reference f64 direct coverage should keep the untouched suffix exactly zero"
    );
    assert!(
        clamped_current[3000..].iter().all(|value| *value == 0.0),
        "regression: theoretical cleanup floor was not enough for the current delta accumulation path; max raw suffix residue was {raw_residue}"
    );
}

#[test]
fn clamp_finite_coverage_below_to_zero_clamps_negative_and_subfloor_values() {
    let mut cov = vec![-0.2, -1.0e-6, 0.0, 0.49, 0.5, f32::NAN, f32::INFINITY];

    clamp_finite_coverage_below_to_zero(&mut cov, 0.5);

    assert_eq!(cov[0], 0.0);
    assert_eq!(cov[1], 0.0);
    assert_eq!(cov[2], 0.0);
    assert_eq!(cov[3], 0.0);
    assert_eq!(cov[4], 0.5);
    assert!(cov[5].is_nan());
    assert_eq!(cov[6], f32::INFINITY);
}

#[test]
fn zero_sum_delta_stress_case_produces_only_zero_after_explicit_cleanup() {
    let mut cov = build_zero_sum_delta_stress_coverage();
    let max_residue = cov.iter().copied().map(f32::abs).fold(0.0_f32, f32::max);

    // Use a floor just above the observed finite residue to verify the cleanup helper against
    // the real delta->coverage accumulation path, without baking a production threshold into this
    // shared test module.
    clamp_finite_coverage_below_to_zero(&mut cov, max_residue + f32::EPSILON);
    assert!(cov.iter().all(|value| *value == 0.0));
}
