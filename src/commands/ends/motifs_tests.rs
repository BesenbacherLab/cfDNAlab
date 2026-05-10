use super::*;
use crate::shared::{
    blacklist::apply_blacklist_mask_to_seq,
    fragment::ends_fragment::{FragmentWithEnds, ResolvedFragmentEnd},
    tiled_run::Tile,
};

fn spec_for_k(k: u8) -> KmerSpec {
    let specs = build_kmer_specs(&[k]).expect("valid k-mer spec");
    specs[&k].clone()
}

fn read_only_motif_context(k_inside: u8) -> TileMotifContext<'static> {
    TileMotifContext {
        reference_start: 0,
        reference_bases: None,
        inside_spec: Some(spec_for_k(k_inside)),
        outside_spec: None,
        inside_codes: None,
        outside_codes: None,
        blacklist_intervals: &[],
        chrom_len: 1_000,
    }
}

fn reference_motif_context(
    seq: &[u8],
    k_inside: Option<u8>,
    k_outside: Option<u8>,
) -> TileMotifContext<'static> {
    reference_motif_context_with_chrom_len(seq, k_inside, k_outside, seq.len() as u64)
}

fn reference_motif_context_with_chrom_len(
    seq: &[u8],
    k_inside: Option<u8>,
    k_outside: Option<u8>,
    chrom_len: u64,
) -> TileMotifContext<'static> {
    let inside_spec = k_inside.map(spec_for_k);
    let outside_spec = k_outside.map(spec_for_k);
    let (inside_codes, outside_codes) = match (inside_spec.as_ref(), outside_spec.as_ref()) {
        (Some(inside_spec), Some(outside_spec)) if inside_spec.k == outside_spec.k => {
            let shared_codes = build_precomputed_reference_codes(Some(inside_spec), seq);
            (shared_codes.clone(), shared_codes)
        }
        _ => (
            build_precomputed_reference_codes(inside_spec.as_ref(), seq),
            build_precomputed_reference_codes(outside_spec.as_ref(), seq),
        ),
    };

    TileMotifContext {
        reference_start: 0,
        reference_bases: Some(seq.to_vec()),
        inside_spec,
        outside_spec,
        inside_codes,
        outside_codes,
        blacklist_intervals: &[],
        chrom_len,
    }
}

fn fragment_with_two_ends(
    left_boundary_pos: u32,
    left_inside: &[u8],
    right_boundary_pos: u32,
    right_inside: &[u8],
) -> FragmentWithEnds {
    FragmentWithEnds {
        tid: 0,
        interval: Interval::new(left_boundary_pos, right_boundary_pos).expect("valid interval"),
        assignment_interval: Interval::new(left_boundary_pos, right_boundary_pos)
            .expect("valid assignment interval"),
        gc_tag: Default::default(),
        left_end: Some(ResolvedFragmentEnd {
            boundary_pos: left_boundary_pos,
            inside_bases: left_inside.to_vec(),
            inside_reference_validation_bp: left_inside.len(),
        }),
        right_end: Some(ResolvedFragmentEnd {
            boundary_pos: right_boundary_pos,
            inside_bases: right_inside.to_vec(),
            inside_reference_validation_bp: right_inside.len(),
        }),
    }
}

#[test]
fn count_fragment_in_window_zero_weight_does_not_create_sparse_window_or_flags() {
    // Arrange: the fragment has two valid read-backed motifs and the window assignment would
    // otherwise count both ends.
    let fragment = fragment_with_two_ends(10, b"AC", 20, b"GT");
    let motif_context = read_only_motif_context(2);
    let mut counts_by_window = EndCountsByWindow::default();

    // Act
    let counted = count_fragment_in_window(
        &mut counts_by_window,
        3,
        Interval::new(10_u64, 20_u64).expect("valid window"),
        &fragment,
        0.0,
        &motif_context,
        KmerSource::Read,
        WindowMotifAssigner::Any,
    )
    .expect("counting should work");

    // Assert: zero-weight motifs should not produce sparse payload entries or counted-end stats.
    assert_eq!(counted, CountedEndFlags::default());
    assert!(counts_by_window.is_empty());
}

#[test]
fn count_fragment_in_window_endpoint_counts_only_left_end_when_window_hits_left_terminal_base() {
    // Arrange: left boundary is at genomic position 10 and right boundary is at 20 for the
    // half-open fragment interval [10, 20).
    //
    // Mental derivation:
    // - endpoint mode counts an end only if its own terminal base lies in the current window
    // - window [10, 12) contains the left terminal base at 10
    // - the right terminal base is `20 - 1 = 19`, which is outside [10, 12)
    // So only the left key should appear, with the provided weight 2.0.
    let fragment = fragment_with_two_ends(10, b"AC", 20, b"GT");
    let motif_context = read_only_motif_context(2);
    let mut counts_by_window = EndCountsByWindow::default();
    let left_key = EncodedEndMotifKey {
        inside_code: motif_context
            .inside_spec
            .as_ref()
            .expect("inside spec")
            .encode_kmer_bytes(b"AC"),
        outside_code: 0,
        reverse_on_decode: false,
    };
    let right_key = EncodedEndMotifKey {
        inside_code: motif_context
            .inside_spec
            .as_ref()
            .expect("inside spec")
            .encode_kmer_bytes(b"GT"),
        outside_code: 0,
        reverse_on_decode: true,
    };

    // Act
    let counted = count_fragment_in_window(
        &mut counts_by_window,
        3,
        Interval::new(10_u64, 12_u64).expect("valid window"),
        &fragment,
        2.0,
        &motif_context,
        KmerSource::Read,
        WindowMotifAssigner::Endpoint,
    )
    .expect("counting should work");

    // Assert
    assert_eq!(
        counted,
        CountedEndFlags {
            left_counted: true,
            right_counted: false,
        }
    );
    let counts = counts_by_window.get(&3).expect("window should be present");
    assert_eq!(counts.counts.get(&left_key), Some(&2.0));
    assert_eq!(counts.counts.get(&right_key), None);
}

#[test]
fn count_fragment_in_window_endpoint_counts_only_right_end_when_window_hits_right_terminal_base() {
    // Arrange:
    // - right boundary is 20, so the right terminal base is `20 - 1 = 19`
    // - endpoint mode should therefore count the right motif in [19, 20), but not the left motif
    let fragment = fragment_with_two_ends(10, b"AC", 20, b"GT");
    let motif_context = read_only_motif_context(2);
    let mut counts_by_window = EndCountsByWindow::default();
    let right_key = EncodedEndMotifKey {
        inside_code: motif_context
            .inside_spec
            .as_ref()
            .expect("inside spec")
            .encode_kmer_bytes(b"GT"),
        outside_code: 0,
        reverse_on_decode: true,
    };

    // Act
    let counted = count_fragment_in_window(
        &mut counts_by_window,
        4,
        Interval::new(19_u64, 20_u64).expect("valid window"),
        &fragment,
        1.25,
        &motif_context,
        KmerSource::Read,
        WindowMotifAssigner::Endpoint,
    )
    .expect("counting should work");

    // Assert
    assert_eq!(
        counted,
        CountedEndFlags {
            left_counted: false,
            right_counted: true,
        }
    );
    let counts = counts_by_window.get(&4).expect("window should be present");
    assert_eq!(counts.counts.len(), 1);
    assert_eq!(counts.counts.get(&right_key), Some(&1.25));
}

#[test]
fn count_fragment_in_window_any_counts_both_ends_in_same_window() {
    // Arrange: `Any` is fragment-centric rather than endpoint-specific.
    //
    // Mental derivation:
    // - once the outer overlap logic has decided this fragment belongs in the window,
    //   `count_fragment_in_window(..., Any)` should count every kept end in that same window
    // - both the left and right motif keys should therefore be present with weight 1.0
    let fragment = fragment_with_two_ends(10, b"AC", 20, b"GT");
    let motif_context = read_only_motif_context(2);
    let mut counts_by_window = EndCountsByWindow::default();
    let left_key = EncodedEndMotifKey {
        inside_code: motif_context
            .inside_spec
            .as_ref()
            .expect("inside spec")
            .encode_kmer_bytes(b"AC"),
        outside_code: 0,
        reverse_on_decode: false,
    };
    let right_key = EncodedEndMotifKey {
        inside_code: motif_context
            .inside_spec
            .as_ref()
            .expect("inside spec")
            .encode_kmer_bytes(b"GT"),
        outside_code: 0,
        reverse_on_decode: true,
    };

    // Act
    let counted = count_fragment_in_window(
        &mut counts_by_window,
        8,
        Interval::new(10_u64, 20_u64).expect("valid window"),
        &fragment,
        1.0,
        &motif_context,
        KmerSource::Read,
        WindowMotifAssigner::Any,
    )
    .expect("counting should work");

    // Assert
    assert_eq!(
        counted,
        CountedEndFlags {
            left_counted: true,
            right_counted: true,
        }
    );
    let counts = counts_by_window.get(&8).expect("window should be present");
    assert_eq!(counts.counts.len(), 2);
    assert_eq!(counts.counts.get(&left_key), Some(&1.0));
    assert_eq!(counts.counts.get(&right_key), Some(&1.0));
}

#[test]
fn encode_inside_code_read_uses_the_resolved_read_bases_directly() {
    // Arrange: in read-backed mode, the inside code should come directly from `inside_bases`
    // rather than from genomic coordinates. So the expected answer is just the codec applied
    // to the literal bytes `AC`.
    let motif_context = read_only_motif_context(2);
    let left_end = ResolvedFragmentEnd {
        boundary_pos: 10,
        inside_bases: b"AC".to_vec(),
        inside_reference_validation_bp: 2,
    };

    // Act
    let code = encode_inside_code(&left_end, EndSide::Left, &motif_context, KmerSource::Read)
        .expect("read-backed inside code should work");

    // Assert
    let spec = motif_context.inside_spec.as_ref().expect("inside spec");
    assert_eq!(code, spec.encode_kmer_bytes(b"AC"));
}

#[test]
fn validate_blacklist_for_read_inside_code_returns_masked_reference_code_for_read_source() {
    // Arrange: the blacklist masks the first inside base at genomic position 2, so the
    // validation code should become the `N` sentinel even though the actual motif would come
    // from the read.
    //
    // Mental derivation:
    // - left inside with `k=2` at boundary 2 reads genomic bases [2, 4)
    // - after masking [2, 3), that span starts with `N`
    // - any k-mer containing `N` must encode as `sentinel_n`
    let inside_spec = spec_for_k(2);
    let mut reference_bases = b"ACGTAC".to_vec();
    let blacklist = [Interval::new(2_u64, 3_u64).expect("valid blacklist")];
    apply_blacklist_mask_to_seq(&mut reference_bases, &blacklist, 0);
    let inside_codes = build_precomputed_reference_codes(Some(&inside_spec), &reference_bases);
    let motif_context = TileMotifContext {
        reference_start: 0,
        reference_bases: Some(reference_bases),
        inside_spec: Some(inside_spec.clone()),
        outside_spec: None,
        inside_codes,
        outside_codes: None,
        blacklist_intervals: &blacklist,
        chrom_len: 6,
    };
    let left_end = ResolvedFragmentEnd {
        boundary_pos: 2,
        inside_bases: b"GT".to_vec(),
        inside_reference_validation_bp: 2,
    };

    // Act
    let code = validate_blacklist_for_read_inside_code(
        &left_end,
        EndSide::Left,
        &inside_spec,
        &motif_context,
    )
    .expect("blacklist validation should work");

    // Assert
    assert_eq!(code, Some(inside_spec.sentinel_n()));
}

#[test]
fn validate_blacklist_for_read_inside_code_ignores_clipped_only_prefix_in_include_at_aligned_boundary_mode() {
    // Arrange: the blacklist masks genomic position 1, but this left end only validates the
    // aligned-overlapping suffix at position 2. The clipped-only prefix must not trigger skipping.
    let inside_spec = spec_for_k(2);
    let mut reference_bases = b"ACGTAC".to_vec();
    let blacklist = [Interval::new(1_u64, 2_u64).expect("valid blacklist")];
    apply_blacklist_mask_to_seq(&mut reference_bases, &blacklist, 0);
    let inside_codes = build_precomputed_reference_codes(Some(&inside_spec), &reference_bases);
    let motif_context = TileMotifContext {
        reference_start: 0,
        reference_bases: Some(reference_bases),
        inside_spec: Some(inside_spec.clone()),
        outside_spec: None,
        inside_codes,
        outside_codes: None,
        blacklist_intervals: &blacklist,
        chrom_len: 6,
    };
    let left_end = ResolvedFragmentEnd {
        boundary_pos: 2,
        inside_bases: b"TG".to_vec(),
        inside_reference_validation_bp: 1,
    };

    // Act
    let code = validate_blacklist_for_read_inside_code(
        &left_end,
        EndSide::Left,
        &inside_spec,
        &motif_context,
    )
    .expect("blacklist validation should work");

    // Assert
    assert_eq!(code, None);
}

#[test]
fn encode_outside_code_left_uses_reference_bases_before_boundary() {
    // Arrange: left outside at boundary 4 with k=2 should read bases [2, 4) = "GT".
    let motif_context = reference_motif_context(b"ACGTAC", None, Some(2));

    // Act
    let code =
        encode_outside_code(4, EndSide::Left, &motif_context).expect("outside code should work");

    // Assert
    let spec = motif_context.outside_spec.as_ref().expect("outside spec");
    assert_eq!(code, spec.encode_kmer_bytes(b"GT"));
}

#[test]
fn encode_outside_code_right_uses_reference_bases_starting_at_boundary() {
    // Arrange: right outside at boundary 2 with k=2 should read bases [2, 4) = "GT".
    let motif_context = reference_motif_context(b"ACGTAC", None, Some(2));

    // Act
    let code =
        encode_outside_code(2, EndSide::Right, &motif_context).expect("outside code should work");

    // Assert
    let spec = motif_context.outside_spec.as_ref().expect("outside spec");
    assert_eq!(code, spec.encode_kmer_bytes(b"GT"));
}

#[test]
fn encode_inside_code_reference_right_uses_reference_bases_ending_at_boundary() {
    // Arrange: right inside at boundary 4 with k=2 should read bases [2, 4) = "GT".
    let motif_context = reference_motif_context(b"ACGTAC", Some(2), None);
    let right_end = ResolvedFragmentEnd {
        boundary_pos: 4,
        inside_bases: vec![],
        inside_reference_validation_bp: 2,
    };

    // Act
    let code = encode_inside_code(
        &right_end,
        EndSide::Right,
        &motif_context,
        KmerSource::Reference,
    )
    .expect("reference-backed inside code should work");

    // Assert
    let spec = motif_context.inside_spec.as_ref().expect("inside spec");
    assert_eq!(code, spec.encode_kmer_bytes(b"GT"));
}

#[test]
fn reference_motif_context_uses_equivalent_codes_for_inside_and_outside_when_k_match() {
    // Arrange: matching `k` values may share an internal table, but the observable contract is
    // simply that the same genomic slice produces the same code for each half.
    let motif_context = reference_motif_context(b"ACGTACGT", Some(2), Some(2));
    let inside_spec = motif_context.inside_spec.as_ref().expect("inside spec");
    let outside_spec = motif_context.outside_spec.as_ref().expect("outside spec");

    // Act
    let inside_code = get_reference_code(
        2,
        inside_spec,
        motif_context.inside_codes.as_deref(),
        &motif_context,
    )
    .expect("inside code");
    let outside_code = get_reference_code(
        2,
        outside_spec,
        motif_context.outside_codes.as_deref(),
        &motif_context,
    )
    .expect("outside code");

    // Assert
    assert_eq!(inside_code, inside_spec.encode_kmer_bytes(b"GT"));
    assert_eq!(outside_code, outside_spec.encode_kmer_bytes(b"GT"));
}

#[test]
fn motif_reference_span_for_tile_extends_full_tile_fetch_by_k_outside_when_aligned() {
    // Mental derivation:
    // - tile fetch is [10, 30)
    // - aligned mode adds no soft-clip expansion
    // - with k_outside=3, preload must cover [10 - 3, 30 + 3) = [7, 33)
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 12, 28, 10, 30).expect("valid tile");

    let reference_span = motif_reference_span_for_tile(&tile, 100, ClipStrategy::Aligned, 5, 3)
        .expect("reference span");

    assert_eq!(
        reference_span,
        Interval::new(7_u64, 33_u64).expect("valid interval")
    );
}

#[test]
fn motif_reference_span_for_tile_keeps_aligned_padding_when_raw_boundary_stays_aligned() {
    // Mental derivation:
    // - tile fetch is [10, 30)
    // - include-at-aligned-boundary adds no soft-clip expansion to the reference preload
    // - with k_outside=3, preload stays [10 - 3, 30 + 3) = [7, 33)
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 12, 28, 10, 30).expect("valid tile");

    let reference_span =
        motif_reference_span_for_tile(&tile, 100, ClipStrategy::IncludeAtAlignedBoundary, 5, 3)
            .expect("reference span");

    assert_eq!(
        reference_span,
        Interval::new(7_u64, 33_u64).expect("valid interval")
    );
}

#[test]
fn motif_reference_span_for_tile_extends_full_tile_fetch_by_k_outside_and_soft_clips_when_raw_boundary_shifts()
 {
    // Mental derivation:
    // - tile fetch is [10, 30)
    // - raw mode adds max_soft_clips on both sides
    // - with k_outside=3 and max_soft_clips=5, pad = 8
    // - preload must therefore cover [10 - 8, 30 + 8) = [2, 38)
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 12, 28, 10, 30).expect("valid tile");

    let reference_span =
        motif_reference_span_for_tile(&tile, 100, ClipStrategy::IncludeAtShiftedBoundary, 5, 3)
            .expect("reference span");

    assert_eq!(
        reference_span,
        Interval::new(2_u64, 38_u64).expect("valid interval")
    );
}

#[test]
fn motif_reference_span_for_tile_clamps_to_chromosome_edges() {
    // Mental derivation:
    // - tile fetch is [4, 18) on a chromosome of length 20
    // - pad is 6 from k_outside=4 plus max_soft_clips=2 in raw mode
    // - unclamped span would be [-2, 24), which must clamp to [0, 20)
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 6, 16, 4, 18).expect("valid tile");

    let reference_span =
        motif_reference_span_for_tile(&tile, 20, ClipStrategy::IncludeAtShiftedBoundary, 2, 4)
            .expect("reference span");

    assert_eq!(
        reference_span,
        Interval::new(0_u64, 20_u64).expect("valid interval")
    );
}

#[test]
fn get_reference_code_errors_when_lookup_escapes_preloaded_span() {
    // Arrange: the preloaded reference slice is [0, 4), but the chromosome extends to 10.
    // Requesting start 3 with k=2 therefore needs [3, 5), which stays within the chromosome
    // yet still escapes the loaded tile slice and must error.
    let motif_context = reference_motif_context_with_chrom_len(b"ACGT", Some(2), None, 10);
    let inside_spec = motif_context.inside_spec.as_ref().expect("inside spec");

    // Act
    let err = get_reference_code(
        3,
        inside_spec,
        motif_context.inside_codes.as_deref(),
        &motif_context,
    )
    .expect_err("lookup outside the preloaded span should fail loudly");

    // Assert
    assert!(
        err.to_string()
            .contains("motif reference lookup escaped preloaded tile span")
    );
}
