use super::*;
use crate::shared::kmers::kmer_codec::build_kmer_specs;

fn spec_for_k(k: u8) -> KmerSpec {
    let specs = build_kmer_specs(&[k]).expect("valid k-mer spec");
    specs[&k].clone()
}

#[test]
fn decode_full_motif_keeps_left_end_in_storage_order() {
    // Arrange: left-end storage order is outside || inside.
    let inside_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);
    let key = EncodedEndMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"GT"),
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: false,
    };

    // Act
    let motif = decode_full_motif(key, Some(&inside_spec), Some(&outside_spec));

    // Assert
    assert_eq!(motif, "ACGT");
}

#[test]
fn decode_full_motif_reverse_complements_right_end() {
    // Arrange: right-end storage order is inside || outside, then reverse-complemented.
    // "AAAC" reverse-complements to "GTTT".
    let inside_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);
    let key = EncodedEndMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"AA"),
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: true,
    };

    // Act
    let motif = decode_full_motif(key, Some(&inside_spec), Some(&outside_spec));

    // Assert
    assert_eq!(motif, "GTTT");
}

#[test]
fn decode_end_motif_counts_collapses_same_orientation_complements_when_requested() {
    // Arrange: decode has already fixed orientation, so collapse compares against the
    // same-orientation complement:
    // - complement("GT") = "CA"
    // - complement("CA") = "GT"
    // - canonical is therefore "CA" for both motifs
    let inside_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"GT"),
            outside_code: 0,
            reverse_on_decode: false,
        },
        1.5,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"CA"),
            outside_code: 0,
            reverse_on_decode: false,
        },
        2.0,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&inside_spec), None, true);

    // Assert
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded.get("_CA"), Some(&3.5));
}

#[test]
fn decode_end_motif_counts_collapses_same_orientation_complements_after_right_end_decode() {
    // Arrange: right-end keys are reverse-complemented during decode first, then collapse is
    // applied on the already oriented motif.
    //
    // Mental derivation:
    // - right-end storage `_AC` decodes to `_GT`
    // - right-end storage `_TG` decodes to `_CA`
    // - `GT` and `CA` are same-orientation complements
    // - both canonicalize to `_CA`
    let inside_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"AC"),
            outside_code: 0,
            reverse_on_decode: true,
        },
        1.25,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"TG"),
            outside_code: 0,
            reverse_on_decode: true,
        },
        2.25,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&inside_spec), None, true);

    // Assert
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded.get("_CA"), Some(&3.5));
}

#[test]
fn decode_end_motif_counts_preserves_outside_inside_order_when_collapsing_combined_motifs() {
    // Arrange: `outside || inside` is the public contract, so collapse must not swap the halves.
    //
    // Mental derivation:
    // - left-end storage outside="G", inside="TA" decodes to "GTA"
    // - right-end storage inside="AT", outside="G" decodes via RC("ATG") to "CAT"
    // - complement("GTA") = "CAT", so both motifs canonicalize to "CAT"
    // - with k_outside=1 and k_inside=2, the final label must stay "C_AT"
    let inside_spec = spec_for_k(2);
    let outside_spec = spec_for_k(1);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"TA"),
            outside_code: outside_spec.encode_kmer_bytes(b"G"),
            reverse_on_decode: false,
        },
        1.25,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"AT"),
            outside_code: outside_spec.encode_kmer_bytes(b"G"),
            reverse_on_decode: true,
        },
        2.25,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&inside_spec), Some(&outside_spec), true);

    // Assert
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded.get("C_AT"), Some(&3.5));
}

#[test]
fn decode_end_motif_counts_preserves_outside_inside_order_for_multi_base_2_plus_2_motifs() {
    // Arrange: this is the same contract as the 1+2 case, but with a wider split.
    //
    // Mental derivation:
    // - left storage outside="GT", inside="AC" decodes to "GTAC"
    // - right storage inside="CA", outside="TG" decodes via RC("CATG") to "CATG"
    // - those are same-orientation complements, so both canonicalize to "CATG" -> "CA_TG"
    //
    // And independently:
    // - left storage outside="TG", inside="CA" decodes to "TGCA"
    // - right storage inside="AC", outside="GT" decodes via RC("ACGT") to "ACGT"
    // - those canonicalize to "ACGT" -> "AC_GT"
    let inside_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"AC"),
            outside_code: outside_spec.encode_kmer_bytes(b"GT"),
            reverse_on_decode: false,
        },
        1.25,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"CA"),
            outside_code: outside_spec.encode_kmer_bytes(b"TG"),
            reverse_on_decode: true,
        },
        2.25,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"CA"),
            outside_code: outside_spec.encode_kmer_bytes(b"TG"),
            reverse_on_decode: false,
        },
        1.0,
    );
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"AC"),
            outside_code: outside_spec.encode_kmer_bytes(b"GT"),
            reverse_on_decode: true,
        },
        1.5,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&inside_spec), Some(&outside_spec), true);

    // Assert
    assert_eq!(decoded.len(), 2);
    assert_eq!(decoded.get("CA_TG"), Some(&3.5));
    assert_eq!(decoded.get("AC_GT"), Some(&2.5));
}

#[test]
fn decode_end_motif_counts_formats_outside_only_labels() {
    // Arrange: with `k_inside = 0`, the public label should still keep the trailing separator.
    let outside_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: 0,
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
    let inside_spec = spec_for_k(2);
    let mut counts = EndMotifCounts::new();
    counts.incr_weighted(
        EncodedEndMotifKey {
            inside_code: inside_spec.encode_kmer_bytes(b"AN"),
            outside_code: 0,
            reverse_on_decode: false,
        },
        1.0,
    );

    // Act
    let decoded = decode_end_motif_counts(&counts, Some(&inside_spec), None, false);

    // Assert
    assert!(decoded.is_empty());
}

#[test]
fn format_end_motif_label_formats_full_motif_as_outside_inside() {
    // Arrange / Act
    let inside_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);
    let label = format_end_motif_label("ACGT", Some(&inside_spec), Some(&outside_spec));

    // Assert
    assert_eq!(label, "AC_GT");
}
