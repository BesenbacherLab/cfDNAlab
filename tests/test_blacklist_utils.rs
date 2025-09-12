#[cfg(test)]
mod tests_merge_intervals {
    use cfdnalab::utils::blacklist::load::merge_intervals;

    #[test]
    fn empty_input() {
        let ivs: Vec<(u64, u64)> = vec![];
        assert!(merge_intervals(ivs).is_empty());
    }

    #[test]
    fn single_interval() {
        let ivs = vec![(100, 200)];
        assert_eq!(merge_intervals(ivs), vec![(100, 200)]);
    }

    #[test]
    fn already_disjoint() {
        let ivs = vec![(10, 20), (30, 40), (50, 60)];
        assert_eq!(
            merge_intervals(ivs.clone()),
            ivs, // should stay exactly the same
        );
    }

    #[test]
    fn overlapping_intervals() {
        // (10, 25) and (20, 40) overlap; (50, 55) is separate
        let ivs = vec![(10, 25), (20, 40), (50, 55)];
        assert_eq!(merge_intervals(ivs), vec![(10, 40), (50, 55)],);
    }

    #[test]
    fn touching_intervals() {
        // Adjacent intervals (end == start) must be coalesced
        let ivs = vec![(0, 10), (10, 20), (20, 30)];
        assert_eq!(merge_intervals(ivs), vec![(0, 30)],);
    }

    #[test]
    fn chain_of_overlaps() {
        // A -> B -> C where each overlaps the next
        let ivs = vec![(1, 5), (4, 8), (7, 12)];
        assert_eq!(merge_intervals(ivs), vec![(1, 12)],);
    }

    #[test]
    fn mixed_sizes_and_overlaps() {
        // Mix of single-base and larger intervals, some overlapping/touching
        //
        // Layout (sorted by start):
        //   (5,6)          – single-base, isolated
        //   (10,100)       – large block
        //   (100,101)      – touches previous -> should merge with (10,100)
        //   (150,160)      – large block
        //   (155,156)      – inside previous -> should merge into (150,160)
        //   (200,201)      – single-base, isolated
        let ivs = vec![
            (5, 6),
            (10, 100),
            (100, 101),
            (150, 160),
            (155, 156),
            (200, 201),
        ];

        assert_eq!(
            merge_intervals(ivs),
            vec![(5, 6), (10, 101), (150, 160), (200, 201)],
        );
    }
}

#[cfg(test)]
mod tests_seq_blacklisting {
    use cfdnalab::utils::blacklist::{apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE};

    #[test]
    fn mask_simple() {
        let mut seq = b"ACGTACGT".to_vec();
        let ivs = vec![(2, 4), (6, 8)]; // mask "GT" and last "GT"
        apply_blacklist_mask_to_seq(&mut seq, &ivs);
        assert_eq!(seq, b"ACXXACXX");
    }

    #[test]
    fn mask_past_end_is_safe() {
        let mut seq = b"AAAA".to_vec();
        let ivs = vec![(2, 10)]; // interval overhangs chromosome
        apply_blacklist_mask_to_seq(&mut seq, &ivs);
        assert_eq!(seq, b"AAXX");
    }

    #[test]
    fn no_intervals_no_change() {
        let original = b"TGCA".to_vec();
        let mut seq = original.clone();
        apply_blacklist_mask_to_seq(&mut seq, &[]);
        assert_eq!(seq, original);
    }

    #[test]
    fn uses_correct_byte() {
        let mut seq = b"GGGG".to_vec();
        apply_blacklist_mask_to_seq(&mut seq, &[(0, 4)]);
        assert!(seq.iter().all(|&b| b == BLACKLIST_BYTE));
    }
}
