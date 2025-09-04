#[cfg(test)]
mod parse_md_tag_tests {
    use cfdna_utils::cfdna_utils::read::parse_md_tag;

    fn check(md: &str, offset: u32, exp_starts: Vec<u32>, exp_ends: Vec<u32>, exp_total: u32) {
        let (starts, ends) = parse_md_tag(md, offset);
        assert_eq!(starts, exp_starts, "Unexpected starts for '{}'", md);
        assert_eq!(ends, exp_ends, "Unexpected ends for '{}'", md);

        // Calculate total number of mismatches across all runs
        let total: u32 = starts.iter().zip(&ends).map(|(&s, &e)| e - s).sum();
        assert_eq!(total, exp_total, "Unexpected total for '{}'", md);
    }

    #[test]
    fn test_only_numbers() {
        // No mismatches
        check("100", 0, vec![], vec![], 0);
        check("0", 0, vec![], vec![], 0);
    }

    #[test]
    fn test_single_mismatch() {
        // One mismatch at pos 10
        check("10A5", 0, vec![10], vec![11], 1);
    }

    #[test]
    fn test_consecutive_mismatches() {
        // Four mismatches in a row starting at pos 5
        check("5AGCT3", 0, vec![5], vec![9], 4);
    }

    #[test]
    fn test_multiple_runs() {
        // Two separate runs: A at pos 3, G at pos 6
        check("3A2G4", 0, vec![3, 6], vec![4, 7], 2);
    }

    #[test]
    fn test_deletion_only() {
        // Deletion of AC at pos 8, no mismatches
        check("8^AC2", 0, vec![], vec![], 0);
    }

    #[test]
    fn test_deletion_and_mismatch() {
        // Deletion at pos 8, then 2 positions, then T mismatch at pos 12
        check("8^AC2T3", 0, vec![12], vec![13], 1);
    }

    #[test]
    fn test_end_mismatch_run() {
        // Mismatch at last base
        check("2T", 0, vec![2], vec![3], 1);
    }

    #[test]
    fn test_empty_string() {
        // Empty MD tag
        check("", 0, vec![], vec![], 0);
    }

    #[test]
    fn test_two_isolated_mismatches() {
        // Two mismatches at pos 2 and 5
        check("2A2C2", 0, vec![2, 5], vec![3, 6], 2);
    }

    #[test]
    fn test_two_consecutive_mismatches() {
        // Consecutive mismatches A and T starting at pos 5
        check("5A0T3", 0, vec![5], vec![7], 2);
    }

    #[test]
    fn test_leading_mismatch() {
        // Leading mismatch at pos 0
        check("0G4", 0, vec![0], vec![1], 1);
    }

    #[test]
    fn test_three_run_then_one() {
        // Run of three mismatches at 3-6, then one at 7
        check("3A0T0C1G2", 0, vec![3, 7], vec![6, 8], 4);
    }
}
