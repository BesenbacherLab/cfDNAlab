#[cfg(test)]
mod tests_seq_blacklisting {
    use crate::shared::{
        blacklist::{apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE},
        interval::Interval,
    };

    fn intervals(entries: &[(u64, u64)]) -> Vec<Interval<u64>> {
        Interval::from_tuples(entries).expect("test intervals should be valid")
    }

    #[test]
    fn mask_simple() {
        let mut seq = b"ACGTACGT".to_vec();
        let ivs = intervals(&[(2, 4), (6, 8)]); // mask "GT" and last "GT"
        apply_blacklist_mask_to_seq(&mut seq, &ivs, 0);
        assert_eq!(seq, b"ACXXACXX");
    }

    #[test]
    fn mask_past_end_is_safe() {
        let mut seq = b"AAAA".to_vec();
        let ivs = intervals(&[(2, 10)]); // interval overhangs chromosome
        apply_blacklist_mask_to_seq(&mut seq, &ivs, 0);
        assert_eq!(seq, b"AAXX");
    }

    #[test]
    fn no_intervals_no_change() {
        let original = b"TGCA".to_vec();
        let mut seq = original.clone();
        apply_blacklist_mask_to_seq(&mut seq, &[], 0);
        assert_eq!(seq, original);
    }

    #[test]
    fn uses_correct_byte() {
        let mut seq = b"GGGG".to_vec();
        let intervals = intervals(&[(0, 4)]);
        apply_blacklist_mask_to_seq(&mut seq, &intervals, 0);
        assert!(seq.iter().all(|&b| b == BLACKLIST_BYTE));
    }

    #[test]
    fn masks_with_offset_slice() {
        let mut seq = b"ACGTACGT".to_vec();
        let ivs = intervals(&[(4, 6)]);
        apply_blacklist_mask_to_seq(&mut seq, &ivs, 2);
        assert_eq!(seq, b"ACXXACGT");
    }
}
