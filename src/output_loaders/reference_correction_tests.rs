use super::*;

// Keep this fixture and its rational expectations identical to the Python and R
// core correction tests. It is the language-parity specification for the math.
const SHARED_COUNTS: [f64; 4] = [2.0, 4.0, 6.0, 8.0];
const SHARED_REFERENCE_FREQUENCIES: [f64; 4] = [1.0 / 8.0, 1.0 / 8.0, 1.0 / 4.0, 1.0 / 2.0];

fn shared_shape() -> TwoSidedMotifShape {
    TwoSidedMotifShape {
        outside_width: 1,
        inside_width: 1,
    }
}

fn shared_parsed_motifs() -> Vec<ParsedEndMotif> {
    ["A_C", "A_G", "T_C", "T_G"]
        .into_iter()
        .map(|label| {
            let (outside, inside) = label.split_once('_').expect("fixture label must split");
            ParsedEndMotif {
                label: label.to_string(),
                outside: outside.to_string(),
                inside: inside.to_string(),
            }
        })
        .collect()
}

fn shared_reference_cache() -> anyhow::Result<SideReferenceCache> {
    let mut cache = SideReferenceCache::new();
    for (label, frequency) in ["AC", "AG", "TC", "TG"]
        .into_iter()
        .zip(SHARED_REFERENCE_FREQUENCIES)
    {
        cache.add_frequency(label, frequency, shared_shape())?;
    }
    cache.finalize_support_counts();
    Ok(cache)
}

fn shared_reference_caches() -> anyhow::Result<BTreeMap<usize, SideReferenceCache>> {
    Ok(BTreeMap::from([(0, shared_reference_cache()?)]))
}

fn shared_counts_matrix() -> anyhow::Result<DenseMatrix<f64>> {
    Ok(DenseMatrix::from_row_major(
        SHARED_COUNTS.to_vec(),
        1,
        SHARED_COUNTS.len(),
    )?)
}

fn assert_values_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (value_index, (&actual_value, &expected_value)) in
        actual.iter().zip(expected).enumerate()
    {
        assert!(
            (actual_value - expected_value).abs() < 1e-12,
            "value {value_index}: expected {expected_value}, got {actual_value}"
        );
    }
}

#[test]
fn joint_core_uses_full_motif_frequencies() -> anyhow::Result<()> {
    let counts = shared_counts_matrix()?;
    let frequencies = DenseMatrix::from_row_major(
        SHARED_REFERENCE_FREQUENCIES.to_vec(),
        1,
        SHARED_REFERENCE_FREQUENCIES.len(),
    )?;
    let motif_labels = shared_parsed_motifs()
        .into_iter()
        .map(|motif| motif.label)
        .collect::<Vec<_>>();

    let corrected = correct_dense_counts(
        &counts,
        &motif_labels,
        &|row_index, motif_index| {
            Ok(frequencies
                .get(row_index, motif_index)
                .copied()
                .unwrap_or(0.0)
                * 4.0)
        },
        UnsupportedReferencePolicy::Error,
    )?;

    // Four positive reference motifs make the uniform frequency 1/4. Relative
    // to that uniform frequency, reference frequencies [1/8, 1/8, 1/4, 1/2]
    // give correction factors [1/2, 1/2, 1, 2] for [AC, AG, TC, TG]. Dividing
    // the original counts [2, 4, 6, 8] by those factors gives [4, 8, 6, 4].
    assert_values_close(corrected.values_row_major(), &[4.0, 8.0, 6.0, 4.0]);
    Ok(())
}

#[test]
fn split_core_multiplies_outside_and_inside_denominators() -> anyhow::Result<()> {
    let counts = shared_counts_matrix()?;
    let parsed_motifs = shared_parsed_motifs();
    let reference_cache = shared_reference_cache()?;
    let motif_labels = parsed_motifs
        .iter()
        .map(|motif| motif.label.clone())
        .collect::<Vec<_>>();

    let corrected = correct_dense_counts(
        &counts,
        &motif_labels,
        &|_, motif_index| Ok(reference_cache.split_denominator(&parsed_motifs[motif_index])),
        UnsupportedReferencePolicy::Error,
    )?;

    // Two positive labels on each side make each side's uniform frequency 1/2.
    // The outside reference frequencies A=1/4 and T=3/4 are therefore 1/2 and
    // 3/2 times uniform. The inside frequencies C=3/8 and G=5/8 are 3/4 and
    // 5/4 times uniform. For [A_C, A_G, T_C, T_G], multiplying the matching
    // side factors gives [3/8, 5/8, 9/8, 15/8]. Dividing original counts
    // [2, 4, 6, 8] by those factors gives [16/3, 32/5, 16/3, 64/15].
    assert_values_close(
        corrected.values_row_major(),
        &[16.0 / 3.0, 32.0 / 5.0, 16.0 / 3.0, 64.0 / 15.0],
    );
    Ok(())
}

#[test]
fn outside_core_aggregates_counts_before_correction() -> anyhow::Result<()> {
    let counts = shared_counts_matrix()?;
    let parsed_motifs = shared_parsed_motifs();
    let side_axis = SideAxisSelection::new(&parsed_motifs, SideMode::Outside, None)?;
    let aggregated = aggregate_dense_side_counts(&counts, &side_axis)?;
    let reference_caches = shared_reference_caches()?;
    let corrected = side_corrected_counts_data(
        &EndMotifCountsData::Dense(aggregated),
        &[0],
        &side_axis.selected_labels,
        SideMode::Outside,
        &reference_caches,
        UnsupportedReferencePolicy::Error,
    )?;
    let EndMotifCountsData::Dense(corrected) = corrected else {
        panic!("dense core fixture must produce dense corrected counts");
    };

    // Counts aggregate to A_=2+4=6 and T_=6+8=14. Two positive outside labels
    // make the uniform frequency 1/2. Relative to uniform, reference
    // frequencies A=1/4 and T=3/4 give factors 1/2 and 3/2. Dividing the
    // aggregated counts by those factors gives [12, 28/3].
    assert_eq!(side_axis.selected_labels, ["A_", "T_"]);
    assert_values_close(corrected.values_row_major(), &[12.0, 28.0 / 3.0]);
    Ok(())
}

#[test]
fn inside_core_aggregates_counts_before_correction() -> anyhow::Result<()> {
    let counts = shared_counts_matrix()?;
    let parsed_motifs = shared_parsed_motifs();
    let side_axis = SideAxisSelection::new(&parsed_motifs, SideMode::Inside, None)?;
    let aggregated = aggregate_dense_side_counts(&counts, &side_axis)?;
    let reference_caches = shared_reference_caches()?;
    let corrected = side_corrected_counts_data(
        &EndMotifCountsData::Dense(aggregated),
        &[0],
        &side_axis.selected_labels,
        SideMode::Inside,
        &reference_caches,
        UnsupportedReferencePolicy::Error,
    )?;
    let EndMotifCountsData::Dense(corrected) = corrected else {
        panic!("dense core fixture must produce dense corrected counts");
    };

    // Counts aggregate to _C=2+6=8 and _G=4+8=12. Two positive inside labels
    // make the uniform frequency 1/2. Relative to uniform, reference
    // frequencies C=3/8 and G=5/8 give factors 3/4 and 5/4. Dividing the
    // aggregated counts by those factors gives [32/3, 48/5].
    assert_eq!(side_axis.selected_labels, ["_C", "_G"]);
    assert_values_close(corrected.values_row_major(), &[32.0 / 3.0, 48.0 / 5.0]);
    Ok(())
}
