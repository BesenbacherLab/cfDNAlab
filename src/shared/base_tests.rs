use super::*;

/// Reference implementation: the original `match` + `to_ascii_uppercase`.
#[inline(always)]
fn encode_base_match(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

#[test]
fn make_canonical_uses_the_lexicographically_smaller_reverse_complement_when_requested() {
    // Arrange:
    // - RC("GT") = "AC"
    // - lexicographically, "AC" < "GT"
    //
    // So canonicalization must return "AC".
    let motif = "GT".to_string();

    // Act
    let canonical = make_canonical(motif, true, false);

    // Assert
    assert_eq!(canonical, "AC");
}

#[test]
fn make_canonical_uses_the_lexicographically_smaller_same_orientation_complement_when_requested() {
    // Arrange:
    // - complement("GTA") = "CAT"
    // - lexicographically, "CAT" < "GTA" because C < G
    //
    // This is the `ends` collapse case: orientation is already fixed, so only the
    // same-orientation complement is compared.
    let motif = "GTA".to_string();

    // Act
    let canonical = make_canonical(motif, false, false);

    // Assert
    assert_eq!(canonical, "CAT");
}

#[test]
fn make_canonical_keeps_a_reverse_complement_palindrome_unchanged() {
    // Arrange / Act / Assert:
    // Palindromic motifs are equal to their reverse complements, so the original string is already
    // canonical.
    assert_eq!(make_canonical("AT".to_string(), true, false), "AT");
}

#[test]
fn lut_equals_match_for_all_bytes() {
    for byte in 0u8..=255 {
        let from_match = encode_base_match(byte);
        let from_lut = encode_base(byte);
        assert_eq!(
            from_match, from_lut,
            "encode_base differs for byte value {byte}"
        );
    }
}
