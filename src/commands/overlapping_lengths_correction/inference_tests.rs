use super::*;
use crate::commands::overlapping_lengths_correction::package::{
    MixtureFitParameters, OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION,
};
use crate::shared::gc_tag::GCTagValue;
use std::collections::BTreeSet;

const COMPLEX_CHUNK_BOUNDARIES: [u32; 4] = [65_536, 131_072, 196_608, 262_144];
const COMPLEX_CONTEXT_END: u32 = 262_600;

struct PlainFragments {
    fragments: VecDeque<FragmentWithSegments>,
}

impl Iterator for PlainFragments {
    type Item = Result<FragmentWithSegments>;

    fn next(&mut self) -> Option<Self::Item> {
        self.fragments.pop_front().map(Ok)
    }
}

fn fragment(start: u32, end: u32) -> Result<FragmentWithSegments> {
    Ok(FragmentWithSegments {
        tid: 0,
        interval: Interval::new(start, end)?,
        segments: None,
        gc_tag: GCTagValue::default(),
        overlap_length_weight: 1.0,
    })
}

fn package() -> Arc<OverlappingLengthsCorrectionPackage> {
    package_with_lookup(&[100.0, 110.0, 121.0], &[1.0, 3.0], 120)
}

/// Build a structurally complete package with a hand-authored lookup for iterator tests.
fn package_with_lookup(
    length_bin_edges: &[f64],
    combined_weights: &[f64],
    maximum_fragment_length: u32,
) -> Arc<OverlappingLengthsCorrectionPackage> {
    assert_eq!(length_bin_edges.len(), combined_weights.len() + 1);
    let length_bin_midpoints = length_bin_edges
        .windows(2)
        .map(|pair| (pair[0] + pair[1]) / 2.0)
        .collect::<Vec<_>>();
    let bin_count = combined_weights.len();
    Arc::new(OverlappingLengthsCorrectionPackage {
        version: OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION,
        length_bin_edges: length_bin_edges.to_vec(),
        length_bin_midpoints,
        length_bin_base_counts: vec![1; bin_count],
        observed_signal_sums: vec![1.0; bin_count],
        raw_depths: vec![1, 2, 3, 4],
        raw_depth_frequencies: vec![1, 1, 1, 1],
        observed_bias: vec![1.0; bin_count],
        first_fitted_bias: vec![1.0; bin_count],
        first_corrected_bias: vec![1.0; bin_count],
        second_fitted_bias: vec![1.0; bin_count],
        target_bias: vec![1.0; bin_count],
        noise_division_factors: vec![1.0; bin_count],
        skew_division_factors: vec![1.0; bin_count],
        mean_shift_division_factors: vec![1.0; bin_count],
        combined_weights: combined_weights.to_vec(),
        initial_fit: MixtureFitParameters {
            scale_multiplier: 1.0,
            skewness: 0.0,
            mean_fragment_length: 110.0,
        },
        refit: MixtureFitParameters {
            scale_multiplier: 1.0,
            skewness: 0.0,
            mean_fragment_length: 110.0,
        },
        minimum_fragment_length: 100,
        maximum_fragment_length,
        minimum_mapq: 30,
        require_proper_pair: false,
        reads_are_fragments: false,
        ignore_gap: false,
        gc_mode: "none".to_string(),
        scaling_enabled: false,
        blacklist_used: false,
    })
}

/// Build several independently varying overlap patterns centered on prefix-chunk boundaries.
///
/// Starts are locally out of order in the same bounded way as paired fragments returned when their
/// second mates are encountered. The segmented fragment omits a gap crossing each boundary, while
/// the other fragments provide changing depth and average fragment length on both sides.
fn complex_boundary_fragments() -> Result<Vec<FragmentWithSegments>> {
    let mut fragments = Vec::new();
    for boundary in COMPLEX_CHUNK_BOUNDARIES {
        fragments.push(fragment(boundary - 140, boundary + 60)?); // 200 bp
        fragments.push(fragment(boundary - 80, boundary + 40)?); // 120 bp

        let mut segmented = fragment(boundary - 110, boundary + 70)?; // 180 bp
        segmented.segments = Some(
            [
                Interval::new(boundary - 110, boundary - 25)?,
                Interval::new(boundary + 15, boundary + 70)?,
            ]
            .into_iter()
            .collect(),
        );
        fragments.push(segmented);

        fragments.push(fragment(boundary - 30, boundary + 130)?); // 160 bp
        fragments.push(fragment(boundary - 50, boundary + 90)?); // 140 bp
    }
    Ok(fragments)
}

/// Run the same fragment stream with a selected test-only prefix storage partition.
fn collect_with_prefix_chunk_size(
    fragments: &[FragmentWithSegments],
    blacklist: &[Interval<u64>],
    prefix_chunk_bases: usize,
) -> Result<Vec<FragmentWithSegments>> {
    let inner = PlainFragments {
        fragments: fragments.iter().cloned().collect(),
    };
    OverlappingLengthWeightIterator::new(
        inner,
        Some(package_with_lookup(
            &[100.0, 130.0, 160.0, 190.0, 221.0],
            &[0.5, 1.25, 2.0, 4.0],
            220,
        )),
        220,
        Interval::new(0, COMPLEX_CONTEXT_END)?,
        blacklist,
    )
    .with_prefix_chunk_bases(prefix_chunk_bases)
    .collect()
}

#[test]
fn fragment_weight_is_base_pair_average_of_positional_lookup_weights() -> Result<()> {
    let inner = PlainFragments {
        fragments: VecDeque::from([
            fragment(100, 200)?,
            fragment(150, 270)?,
            fragment(400, 500)?,
        ]),
    };
    let mut weighted = OverlappingLengthWeightIterator::new(
        inner,
        Some(package()),
        120,
        Interval::new(0, 700)?,
        &[],
    );
    let first = weighted.next().context("first weighted fragment")??;
    assert_eq!(first.interval.as_tuple(), (100, 200));
    assert!((first.overlap_length_weight - 2.0).abs() < 1.0e-12);
    let second = weighted.next().context("second weighted fragment")??;
    assert_eq!(second.interval.as_tuple(), (150, 270));
    assert!((second.overlap_length_weight - 3.0).abs() < 1.0e-12);
    Ok(())
}

#[test]
fn blacklisted_bases_are_excluded_from_fragment_weight_average() -> Result<()> {
    let inner = PlainFragments {
        fragments: VecDeque::from([fragment(100, 200)?, fragment(150, 270)?]),
    };
    let blacklist = [Interval::new(150_u64, 200_u64)?];
    let mut weighted = OverlappingLengthWeightIterator::new(
        inner,
        Some(package()),
        120,
        Interval::new(0, 400)?,
        &blacklist,
    );
    let first = weighted.next().context("first weighted fragment")??;
    assert!((first.overlap_length_weight - 1.0).abs() < 1.0e-12);
    Ok(())
}

#[test]
fn prefix_chunks_do_not_change_weights_at_chunk_boundary() -> Result<()> {
    let inner = PlainFragments {
        fragments: VecDeque::from([fragment(65_500, 65_600)?, fragment(65_530, 65_650)?]),
    };
    let mut weighted = OverlappingLengthWeightIterator::new(
        inner,
        Some(package()),
        120,
        Interval::new(0, 65_800)?,
        &[],
    );
    let first = weighted.next().context("first weighted fragment")??;
    // The first 30 bases have only the 100 bp fragment and weight 1. The remaining 70 bases
    // overlap both fragments, average 110 bp, and have weight 3. The chunk boundary at 65,536 lies
    // inside the second run, so the expected fragment mean is (30 * 1 + 70 * 3) / 100 = 2.4
    assert!((first.overlap_length_weight - 2.4).abs() < 1.0e-12);
    Ok(())
}

#[test]
fn complex_stream_is_bitwise_identical_with_64_and_256_kib_prefix_chunks() -> Result<()> {
    // Arrange
    assert_eq!(PREFIX_CHUNK_BASES, 65_536);
    let fragments = complex_boundary_fragments()?;
    let mut blacklist = Vec::new();
    for (boundary_index, boundary) in COMPLEX_CHUNK_BOUNDARIES.into_iter().enumerate() {
        // Keep three boundaries observable so a missing or duplicated boundary base changes the
        // result. Mask across one boundary to exercise exclusion while chunks are rotated.
        let first_blacklist = if boundary_index == 1 {
            Interval::new(
                u64::from(boundary - 12),
                u64::from(boundary + 18),
            )?
        } else {
            Interval::new(
                u64::from(boundary - 70),
                u64::from(boundary - 55),
            )?
        };
        blacklist.push(first_blacklist);
        blacklist.push(Interval::new(
            u64::from(boundary + 45),
            u64::from(boundary + 58),
        )?);
    }
    for boundary in COMPLEX_CHUNK_BOUNDARIES {
        assert!(
            fragments
                .iter()
                .any(|fragment| fragment.start() < boundary && boundary < fragment.end()),
            "fixture must contain a fragment crossing boundary {boundary}"
        );
    }

    // Act
    let with_64_kib_chunks =
        collect_with_prefix_chunk_size(&fragments, &blacklist, PREFIX_CHUNK_BASES)?;
    let with_256_kib_chunks = collect_with_prefix_chunk_size(&fragments, &blacklist, 262_144)?;

    // Assert
    assert_eq!(with_64_kib_chunks.len(), fragments.len());
    assert_eq!(with_256_kib_chunks.len(), fragments.len());
    let mut distinct_weight_bits = BTreeSet::new();
    for (fragment_index, ((input, with_64_kib), with_256_kib)) in fragments
        .iter()
        .zip(&with_64_kib_chunks)
        .zip(&with_256_kib_chunks)
        .enumerate()
    {
        assert_eq!(
            with_64_kib.interval, input.interval,
            "64 KiB output order changed at fragment {fragment_index}"
        );
        assert_eq!(
            with_256_kib.interval, input.interval,
            "256 KiB output order changed at fragment {fragment_index}"
        );
        assert_eq!(with_64_kib.segments, input.segments);
        assert_eq!(with_256_kib.segments, input.segments);
        assert_eq!(
            with_64_kib.overlap_length_weight.to_bits(),
            with_256_kib.overlap_length_weight.to_bits(),
            "chunk size changed the weight of fragment {fragment_index}"
        );
        distinct_weight_bits.insert(with_64_kib.overlap_length_weight.to_bits());
    }
    assert!(
        distinct_weight_bits.len() >= 4,
        "fixture should exercise several different inferred weights"
    );
    Ok(())
}

#[test]
fn weighted_iterator_preserves_input_order_when_later_fragment_is_safe_first() -> Result<()> {
    // Fragment B has an earlier end than A. Once C raises the largest observed start to 410,
    // the exclusive safe boundary reaches 290, making B safe while A remains unsafe. D advances
    // that boundary to 330 so the FIFO front A can finally be returned before B
    let input_intervals = [(200, 320), (250, 280), (410, 520), (450, 550)];
    let inner = PlainFragments {
        fragments: input_intervals
            .iter()
            .map(|&(start, end)| fragment(start, end))
            .collect::<Result<VecDeque<_>>>()?,
    };
    let weighted = OverlappingLengthWeightIterator::new(
        inner,
        Some(package()),
        120,
        Interval::new(0, 700)?,
        &[],
    );

    let output_intervals = weighted
        .map(|result| result.map(|fragment| fragment.interval.as_tuple()))
        .collect::<Result<Vec<_>>>()?;

    assert_eq!(output_intervals, input_intervals);
    Ok(())
}
