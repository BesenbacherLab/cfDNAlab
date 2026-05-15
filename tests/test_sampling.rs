#[cfg(test)]
mod tests_sample_starts_in_core {
    use anyhow::Result;
    use cfdnalab::shared::sampling::sample_starts_per_chrom;
    use fxhash::FxHashMap;
    use rand::{SeedableRng, rngs::StdRng};

    fn fmap(pairs: &[(&str, usize)]) -> FxHashMap<String, usize> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect::<FxHashMap<_, _>>()
    }

    #[test]
    fn returns_empty_when_no_possible_positions() -> Result<()> {
        let chrom_sizes = fmap(&[("chrTiny", 50), ("chrSmall", 99)]);
        let mut rng = StdRng::seed_from_u64(42);
        let out = sample_starts_per_chrom(&mut rng, &chrom_sizes, 10_000, 100)?;
        assert!(out.is_empty());
        Ok(())
    }

    #[test]
    fn returns_empty_when_n_samples_zero() -> Result<()> {
        let chrom_sizes = fmap(&[("chr1", 1_000)]);
        let mut rng = StdRng::seed_from_u64(1);
        let out = sample_starts_per_chrom(&mut rng, &chrom_sizes, 0, 100)?;
        assert!(out.is_empty());
        Ok(())
    }

    #[test]
    fn exact_apportionment_even_split() -> Result<()> {
        // Two chromosomes with equal possible starts -> equal quotas.
        // L = 101, max_window_len = 2 -> possible = 100 each; total = 200.
        let chrom_sizes = fmap(&[("chrA", 101), ("chrB", 101)]);
        let mut rng = StdRng::seed_from_u64(7);
        let out = sample_starts_per_chrom(&mut rng, &chrom_sizes, 200, 2)?;
        assert_eq!(out["chrA"].len(), 100);
        assert_eq!(out["chrB"].len(), 100);

        // Check bounds and uniqueness+sortedness
        for (name, starts) in &out {
            assert!(
                starts.windows(2).all(|w| w[0] < w[1]),
                "duplicates or not sorted in {name}"
            );
            let possible = 101 - 2 + 1; // 100 valid starts -> indices 0..=99
            assert!(*starts.last().unwrap() < possible);
        }
        Ok(())
    }

    #[test]
    fn hamilton_rounding_prefers_larger_fraction() -> Result<()> {
        // possible: chrBig=100, chrSmall=50; total=150; n=10
        // exact shares: 6.666..., 3.333... -> floor [6,3], remaining=1 -> add to largest remainder -> [7,3]
        let chrom_sizes = fmap(&[("chrBig", 101), ("chrSmall", 51)]); // possible=100,50
        let mut rng = StdRng::seed_from_u64(99);
        let out = sample_starts_per_chrom(&mut rng, &chrom_sizes, 10, 2)?;
        assert_eq!(out["chrBig"].len(), 7);
        assert_eq!(out["chrSmall"].len(), 3);

        // No duplicates + sorted
        assert!(out["chrBig"].windows(2).all(|w| w[0] < w[1]));
        assert!(out["chrSmall"].windows(2).all(|w| w[0] < w[1]));
        Ok(())
    }

    #[test]
    fn skips_chromosomes_shorter_than_max_window_len() -> Result<()> {
        let chrom_sizes = fmap(&[("chr1", 1000), ("chrTooShort", 5)]);
        let mut rng = StdRng::seed_from_u64(123);
        let out = sample_starts_per_chrom(&mut rng, &chrom_sizes, 100, 100)?;
        assert!(out.contains_key("chr1"));
        assert!(!out.contains_key("chrTooShort"));
        Ok(())
    }

    #[test]
    fn total_samples_equal_to_n_samples() -> Result<()> {
        let chrom_sizes = fmap(&[
            ("chr1", 10_001), // possible 10_000
            ("chr2", 5_501),  // possible  5_500
            ("chr3", 1_001),  // possible  1_000
        ]);
        let mut rng = StdRng::seed_from_u64(2024);
        let n = 1234;
        let out = sample_starts_per_chrom(&mut rng, &chrom_sizes, n, 2)?;
        let total: usize = out.values().map(|v| v.len()).sum();
        assert_eq!(total, n);
        Ok(())
    }

    #[test]
    fn deterministic_with_seeded_rng() -> Result<()> {
        let chrom_sizes = fmap(&[("chrX", 1_000_000), ("chrY", 800_000)]);
        let mut rng1 = StdRng::seed_from_u64(777);
        let mut rng2 = StdRng::seed_from_u64(777);

        let a = sample_starts_per_chrom(&mut rng1, &chrom_sizes, 10_000, 1000)?;
        let b = sample_starts_per_chrom(&mut rng2, &chrom_sizes, 10_000, 1000)?;
        assert_eq!(a, b);
        Ok(())
    }

    #[test]
    fn indices_within_valid_range() -> Result<()> {
        let chrom_sizes = fmap(&[("chrZ", 10_000)]);
        let mut rng = StdRng::seed_from_u64(31337);
        let max_window_len = 123;
        let possible = 10_000 - max_window_len + 1; // inclusive end range
        let out = sample_starts_per_chrom(&mut rng, &chrom_sizes, 500, max_window_len)?;
        let starts = &out["chrZ"];
        assert!(starts.iter().all(|&s| s < possible));
        Ok(())
    }
}

#[cfg(test)]
mod tests_tile_samplers {
    use anyhow::Result;
    use cfdnalab::shared::sampling::{sample_starts_in_core, sampling_density};
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
