#[cfg(test)]
mod tests_tile_samplers {
    use crate::shared::sampling::{sample_starts_in_core, sampling_density};
    use anyhow::Result;
    use fxhash::FxHashMap;
    use rand::{SeedableRng, rngs::StdRng};

    fn fmap(pairs: &[(&str, usize)]) -> FxHashMap<String, usize> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect::<FxHashMap<_, _>>()
    }

    #[test]
    fn sampling_density_zero_when_no_valid_starts() -> Result<()> {
        let chrom_sizes = fmap(&[("chrTiny", 50), ("chrSmall", 80)]);
        let density = sampling_density(&chrom_sizes, 100, 1_000);
        assert_eq!(density, 0.0);
        Ok(())
    }

    #[test]
    fn returns_empty_core_when_requested_samples_are_zero() -> Result<()> {
        let chrom_sizes = fmap(&[("chr1", 1_000)]);
        let density = sampling_density(&chrom_sizes, 100, 0);
        assert_eq!(density, 0.0);

        let mut rng = StdRng::seed_from_u64(1);
        let starts = sample_starts_in_core(&mut rng, 0, 100, 1_000, 100, density);
        assert!(starts.is_empty());
        Ok(())
    }

    #[test]
    fn sampling_density_matches_ratio_of_possible_starts() -> Result<()> {
        // Possible starts: chrLong=11, chrMedium=6 -> total 17; n_samples=34 -> density=2.0
        let chrom_sizes = fmap(&[("chrLong", 20), ("chrMedium", 15)]);
        let density = sampling_density(&chrom_sizes, 10, 34);
        assert_eq!(density, 2.0);
        Ok(())
    }

    #[test]
    fn returns_empty_core_when_no_valid_positions() -> Result<()> {
        // Core end < start or fragment longer than chromosome -> nothing to sample
        let mut rng = StdRng::seed_from_u64(42);
        let empty_core = sample_starts_in_core(&mut rng, 50, 40, 1_000, 100, 0.5);
        let too_long = sample_starts_in_core(&mut rng, 0, 10, 5, 10, 0.5);
        assert!(empty_core.is_empty());
        assert!(too_long.is_empty());
        Ok(())
    }

    #[test]
    fn samples_entire_core_when_estimate_exceeds_available() -> Result<()> {
        let mut rng = StdRng::seed_from_u64(5);
        let starts = sample_starts_in_core(&mut rng, 100, 110, 200, 5, 5.0);
        let expected: Vec<usize> = (100..110).collect();
        assert_eq!(starts, expected);
        Ok(())
    }

    #[test]
    fn limits_starts_to_chromosome_end() -> Result<()> {
        let mut rng = StdRng::seed_from_u64(9);
        let starts = sample_starts_in_core(&mut rng, 80, 120, 100, 6, 10.0);
        let expected: Vec<usize> = (80..95).collect();
        assert_eq!(starts, expected);
        Ok(())
    }

    #[test]
    fn core_sampling_is_deterministic_with_seeded_rng() -> Result<()> {
        let mut rng1 = StdRng::seed_from_u64(777);
        let mut rng2 = StdRng::seed_from_u64(777);
        let a = sample_starts_in_core(&mut rng1, 10, 50, 1_000, 12, 0.5);
        let b = sample_starts_in_core(&mut rng2, 10, 50, 1_000, 12, 0.5);
        assert_eq!(a, b);
        Ok(())
    }

    #[test]
    fn core_sampling_obeys_density_less_than_one() -> Result<()> {
        let mut rng = StdRng::seed_from_u64(1);
        let starts = sample_starts_in_core(&mut rng, 0, 100, 1_000, 20, 0.25);
        assert!(starts.len() <= 25);
        Ok(())
    }
}
