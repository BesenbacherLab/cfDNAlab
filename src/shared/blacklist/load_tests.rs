#[cfg(test)]
mod tests_load_blacklists {
    use crate::shared::{blacklist::load::load_blacklists, interval::Interval};
    use anyhow::Result;
    use tempfile::NamedTempFile;

    fn write_bed(lines: &[&str]) -> Result<NamedTempFile> {
        let mut file = NamedTempFile::new()?;
        use std::io::Write;
        for line in lines {
            writeln!(file, "{}", line)?;
        }
        Ok(file)
    }

    #[test]
    fn should_filter_by_min_size_and_whitelist() -> Result<()> {
        // Arrange
        let bed = write_bed(&["chr1\t0\t3", "chr1\t10\t20", "chr2\t5\t30"])?;
        let whitelist = vec!["chr1".to_string()];

        // Act
        let map = load_blacklists(&[bed.path()], 5, 0, Some(whitelist.as_slice()))?;

        // Assert
        assert_eq!(
            map.get("chr1").unwrap().as_slice(),
            Interval::from_tuples(&[(10, 20)])?.as_slice()
        );
        assert!(!map.contains_key("chr2"));
        Ok(())
    }

    #[test]
    fn should_expand_by_halo_before_merging() -> Result<()> {
        // Arrange
        let bed = write_bed(&["chrX\t100\t110", "chrX\t112\t120"])?;

        // Act
        let map = load_blacklists(&[bed.path()], 1, 2, None)?;

        // Assert
        assert_eq!(
            map.get("chrX").unwrap().as_slice(),
            Interval::from_tuples(&[(98, 122)])?.as_slice()
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests_merge_intervals {
    use crate::shared::blacklist::load::merge_intervals;
    use crate::shared::interval::Interval;

    #[test]
    fn empty_input() -> anyhow::Result<()> {
        let ivs: Vec<Interval<u64>> = vec![];
        assert!(merge_intervals(ivs)?.is_empty());
        Ok(())
    }

    #[test]
    fn single_interval() -> anyhow::Result<()> {
        let ivs = vec![Interval::new(100, 200)?];
        assert_eq!(merge_intervals(ivs)?, vec![Interval::new(100, 200)?]);
        Ok(())
    }

    #[test]
    fn already_disjoint() -> anyhow::Result<()> {
        let ivs = vec![
            Interval::new(10, 20)?,
            Interval::new(30, 40)?,
            Interval::new(50, 60)?,
        ];
        assert_eq!(
            merge_intervals(ivs.clone())?,
            ivs, // should stay exactly the same
        );
        Ok(())
    }

    #[test]
    fn overlapping_intervals() -> anyhow::Result<()> {
        // (10, 25) and (20, 40) overlap; (50, 55) is separate
        let ivs = vec![
            Interval::new(10, 25)?,
            Interval::new(20, 40)?,
            Interval::new(50, 55)?,
        ];
        assert_eq!(
            merge_intervals(ivs)?,
            vec![Interval::new(10, 40)?, Interval::new(50, 55)?],
        );
        Ok(())
    }

    #[test]
    fn touching_intervals() -> anyhow::Result<()> {
        // Adjacent intervals (end == start) must be coalesced
        let ivs = vec![
            Interval::new(0, 10)?,
            Interval::new(10, 20)?,
            Interval::new(20, 30)?,
        ];
        assert_eq!(merge_intervals(ivs)?, vec![Interval::new(0, 30)?]);
        Ok(())
    }

    #[test]
    fn chain_of_overlaps() -> anyhow::Result<()> {
        // A -> B -> C where each overlaps the next
        let ivs = vec![
            Interval::new(1, 5)?,
            Interval::new(4, 8)?,
            Interval::new(7, 12)?,
        ];
        assert_eq!(merge_intervals(ivs)?, vec![Interval::new(1, 12)?]);
        Ok(())
    }

    #[test]
    fn mixed_sizes_and_overlaps() -> anyhow::Result<()> {
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
            Interval::new(5, 6)?,
            Interval::new(10, 100)?,
            Interval::new(100, 101)?,
            Interval::new(150, 160)?,
            Interval::new(155, 156)?,
            Interval::new(200, 201)?,
        ];

        assert_eq!(
            merge_intervals(ivs)?,
            vec![
                Interval::new(5, 6)?,
                Interval::new(10, 101)?,
                Interval::new(150, 160)?,
                Interval::new(200, 201)?,
            ],
        );
        Ok(())
    }
}
