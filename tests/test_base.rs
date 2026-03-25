#[cfg(test)]
mod tests {
    use cfdnalab::shared::base::*;

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
    fn lut_equals_match_for_all_bytes() {
        // Human verification status: unverified
        for byte in 0u8..=255 {
            let from_match = encode_base_match(byte);
            let from_lut = encode_base(byte);
            assert_eq!(
                from_match, from_lut,
                "encode_base differs for byte value {byte}"
            );
        }
    }
}
