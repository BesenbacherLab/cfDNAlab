use super::*;

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
fn make_canonical_uses_the_lexicographically_smaller_same_orientation_complement_when_requested()
{
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
