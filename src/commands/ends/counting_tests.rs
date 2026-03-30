use super::*;
use crate::shared::kmers::kmer_codec::build_kmer_specs;

fn spec_for_k(k: u8) -> KmerSpec {
    let specs = build_kmer_specs(&[k]).expect("valid k-mer spec");
    specs[&k].clone()
}

#[test]
fn decode_full_motif_keeps_left_end_in_storage_order() {
    // Arrange: left-end storage order is outside || within.
    let within_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);
    let key = EncodedEndMotifKey {
        within_code: within_spec.encode_kmer_bytes(b"GT"),
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: false,
    };

    // Act
    let motif = decode_full_motif(key, Some(&within_spec), Some(&outside_spec));

    // Assert
    assert_eq!(motif, "ACGT");
}

#[test]
fn decode_full_motif_reverse_complements_right_end() {
    // Arrange: right-end storage order is within || outside, then reverse-complemented.
    // "AAAC" reverse-complements to "GTTT".
    let within_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);
    let key = EncodedEndMotifKey {
        within_code: within_spec.encode_kmer_bytes(b"AA"),
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: true,
    };

    // Act
    let motif = decode_full_motif(key, Some(&within_spec), Some(&outside_spec));

    // Assert
    assert_eq!(motif, "GTTT");
}

#[test]
fn decode_end_motif_counts_collapses_reverse_complements_when_requested() {
    // Arrange: "GT" and "AC" are reverse complements and both canonicalize to "AC".
    let within_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            within_code: within_spec.encode_kmer_bytes(b"GT"),
            outside_code: 0,
            reverse_on_decode: false,
        },
        1.5,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            within_code: within_spec.encode_kmer_bytes(b"AC"),
            outside_code: 0,
            reverse_on_decode: false,
        },
        2.0,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&within_spec), None, true);

    // Assert
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded.get("_AC"), Some(&3.5));
}

#[test]
fn decode_end_motif_counts_collapses_right_end_reverse_complements_when_requested() {
    // Arrange: for right-end keys, decode first reverse-complements the stored within bases.
    //
    // Mental derivation:
    // - right-end storage `_GT` decodes to `_AC`
    // - right-end storage `_AC` decodes to `_GT`
    // - with complement collapsing enabled, both labels canonicalize to `_AC`
    let within_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            within_code: within_spec.encode_kmer_bytes(b"GT"),
            outside_code: 0,
            reverse_on_decode: true,
        },
        1.25,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            within_code: within_spec.encode_kmer_bytes(b"AC"),
            outside_code: 0,
            reverse_on_decode: true,
        },
        2.25,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&within_spec), None, true);

    // Assert
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded.get("_AC"), Some(&3.5));
}

#[test]
fn decode_end_motif_counts_formats_outside_only_labels() {
    // Arrange: with `k_within = 0`, the public label should still keep the trailing separator.
    let outside_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            within_code: 0,
            outside_code: outside_spec.encode_kmer_bytes(b"AC"),
            reverse_on_decode: false,
        },
        1.25,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, None, Some(&outside_spec), false);

    // Assert
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded.get("AC_"), Some(&1.25));
}

#[test]
fn decode_end_motif_counts_drops_motifs_with_n() {
    // Arrange: sentinel-N decodes to an N-containing motif and should be dropped.
    let within_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            within_code: within_spec.encode_kmer_bytes(b"AN"),
            outside_code: 0,
            reverse_on_decode: false,
        },
        1.0,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&within_spec), None, false);

    // Assert
    assert!(decoded.is_empty());
}

#[test]
fn format_end_motif_label_formats_full_motif_as_outside_within() {
    // Arrange / Act
    let within_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);
    let label = format_end_motif_label("ACGT", Some(&within_spec), Some(&outside_spec));

    // Assert
    assert_eq!(label, "AC_GT");
}
