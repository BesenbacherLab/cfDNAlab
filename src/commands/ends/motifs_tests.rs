use super::*;
use crate::shared::fragment::ends_fragment::{FragmentWithEnds, ResolvedFragmentEnd};
use std::sync::Arc;

fn spec_for_k(k: u8) -> KmerSpec {
    let specs = build_kmer_specs(&[k]).expect("valid k-mer spec");
    specs[&k].clone()
}

fn read_only_motif_context(k_within: u8) -> TileMotifContext<'static> {
    TileMotifContext {
        chrom: "chr1".to_string(),
        ref_2bit: None,
        fetch_start: 0,
        reference_bases: None,
        within_spec: Some(spec_for_k(k_within)),
        outside_spec: None,
        within_codes: None,
        outside_codes: None,
        blacklist_intervals: &[],
        chrom_len: 1_000,
    }
}

fn reference_motif_context(
    seq: &[u8],
    k_within: Option<u8>,
    k_outside: Option<u8>,
) -> TileMotifContext<'static> {
    let within_spec = k_within.map(spec_for_k);
    let outside_spec = k_outside.map(spec_for_k);
    let (within_codes, outside_codes) = match (within_spec.as_ref(), outside_spec.as_ref()) {
        (Some(within_spec), Some(outside_spec)) if within_spec.k == outside_spec.k => {
            let shared_codes = build_precomputed_reference_codes(Some(within_spec), seq);
            (shared_codes.clone(), shared_codes)
        }
        _ => (
            build_precomputed_reference_codes(within_spec.as_ref(), seq),
            build_precomputed_reference_codes(outside_spec.as_ref(), seq),
        ),
    };

    TileMotifContext {
        chrom: "chr1".to_string(),
        ref_2bit: None,
        fetch_start: 0,
        reference_bases: Some(seq.to_vec()),
        within_spec,
        outside_spec,
        within_codes,
        outside_codes,
        blacklist_intervals: &[],
        chrom_len: seq.len() as u64,
    }
}

fn fragment_with_two_ends(
    left_boundary_pos: u32,
    left_within: &[u8],
    right_boundary_pos: u32,
    right_within: &[u8],
) -> FragmentWithEnds {
    FragmentWithEnds {
        tid: 0,
        interval: Interval::new(left_boundary_pos, right_boundary_pos).expect("valid interval"),
        assignment_interval: Interval::new(left_boundary_pos, right_boundary_pos)
            .expect("valid assignment interval"),
        gc_tag: Default::default(),
        left_end: Some(ResolvedFragmentEnd {
            boundary_pos: left_boundary_pos,
            within_bases: left_within.to_vec(),
        }),
        right_end: Some(ResolvedFragmentEnd {
            boundary_pos: right_boundary_pos,
            within_bases: right_within.to_vec(),
        }),
    }
}

#[test]
fn count_fragment_in_window_endpoint_counts_only_left_end_when_window_hits_left_terminal_base() {
    // Arrange
    let fragment = fragment_with_two_ends(10, b"AC", 20, b"GT");
    let motif_context = read_only_motif_context(2);
    let mut counts_by_window = EndCountsByWindow::default();
    let left_key = EncodedEndMotifKey {
        within_code: motif_context
            .within_spec
            .as_ref()
            .expect("within spec")
            .encode_kmer_bytes(b"AC"),
        outside_code: 0,
        reverse_on_decode: false,
    };
    let right_key = EncodedEndMotifKey {
        within_code: motif_context
            .within_spec
            .as_ref()
            .expect("within spec")
            .encode_kmer_bytes(b"GT"),
        outside_code: 0,
        reverse_on_decode: true,
    };

    // Act
    count_fragment_in_window(
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
    let counts = counts_by_window.get(&3).expect("window should be present");
    assert_eq!(counts.counts.get(&left_key), Some(&2.0));
    assert_eq!(counts.counts.get(&right_key), None);
}

#[test]
fn count_fragment_in_window_endpoint_counts_only_right_end_when_window_hits_right_terminal_base() {
    // Arrange
    let fragment = fragment_with_two_ends(10, b"AC", 20, b"GT");
    let motif_context = read_only_motif_context(2);
    let mut counts_by_window = EndCountsByWindow::default();
    let right_key = EncodedEndMotifKey {
        within_code: motif_context
            .within_spec
            .as_ref()
            .expect("within spec")
            .encode_kmer_bytes(b"GT"),
        outside_code: 0,
        reverse_on_decode: true,
    };

    // Act
    count_fragment_in_window(
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
    let counts = counts_by_window.get(&4).expect("window should be present");
    assert_eq!(counts.counts.len(), 1);
    assert_eq!(counts.counts.get(&right_key), Some(&1.25));
}

#[test]
fn count_fragment_in_window_any_counts_both_ends_in_same_window() {
    // Arrange
    let fragment = fragment_with_two_ends(10, b"AC", 20, b"GT");
    let motif_context = read_only_motif_context(2);
    let mut counts_by_window = EndCountsByWindow::default();

    // Act
    count_fragment_in_window(
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
    let counts = counts_by_window.get(&8).expect("window should be present");
    assert_eq!(counts.counts.len(), 2);
}

#[test]
fn encode_outside_code_left_uses_reference_bases_before_boundary() {
    // Arrange: left outside at boundary 4 with k=2 should read bases [2, 4) = "GT".
    let motif_context = reference_motif_context(b"ACGTAC", None, Some(2));

    // Act
    let code = encode_outside_code(4, EndSide::Left, &motif_context)
        .expect("outside code should work");

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
fn encode_within_code_reference_right_uses_reference_bases_ending_at_boundary() {
    // Arrange: right within at boundary 4 with k=2 should read bases [2, 4) = "GT".
    let motif_context = reference_motif_context(b"ACGTAC", Some(2), None);
    let right_end = ResolvedFragmentEnd {
        boundary_pos: 4,
        within_bases: vec![],
    };

    // Act
    let code = encode_within_code(
        &right_end,
        EndSide::Right,
        &motif_context,
        KmerSource::Reference,
    )
    .expect("reference-backed within code should work");

    // Assert
    let spec = motif_context.within_spec.as_ref().expect("within spec");
    assert_eq!(code, spec.encode_kmer_bytes(b"GT"));
}

#[test]
fn reference_motif_context_shares_code_table_when_within_and_outside_k_match() {
    // Arrange
    let motif_context = reference_motif_context(b"ACGTACGT", Some(2), Some(2));

    // Act
    let within_codes = motif_context.within_codes.as_ref().expect("within codes");
    let outside_codes = motif_context.outside_codes.as_ref().expect("outside codes");

    // Assert
    assert!(Arc::ptr_eq(within_codes, outside_codes));
}
