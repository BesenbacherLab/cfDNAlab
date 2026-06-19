use super::*;

fn chromosomes() -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
        chromosomes_file: None,
    }
}

#[test]
fn ref_gc_bias_config_new_uses_cli_defaults() {
    // Arrange
    let ref_2bit = PathBuf::from("ref.2bit");
    let output_dir = PathBuf::from("out");
    let chromosomes = chromosomes();

    // Act
    let config = RefGCBiasConfig::new(ref_2bit.clone(), output_dir.clone(), chromosomes.clone());

    // Assert
    assert_eq!(config.ref_genome.ref_2bit, ref_2bit);
    assert_eq!(config.output_dir, output_dir);
    assert_eq!(config.chromosomes, chromosomes);
    assert_eq!(config.output_prefix, "");
    assert_eq!(
        config.n_threads,
        crate::shared::thread_pool::default_thread_count()
    );
    assert_eq!(config.n_positions, DEFAULT_N_POSITIONS);
    assert_eq!(config.seed, None);
    assert_eq!(config.windows, RefGCWindowsArgs::default());
    assert_eq!(config.blacklist, None);
    assert_eq!(config.fragment_lengths, FragmentLengthArgs::default());
    assert_eq!(config.end_offset, DEFAULT_END_OFFSET);
    assert!(!config.skip_interpolation);
    assert_eq!(config.smoothing_sigma, DEFAULT_SMOOTHING_SIGMA);
    assert_eq!(config.smoothing_radius, DEFAULT_SMOOTHING_RADIUS);
    assert!(!config.skip_smoothing);
    assert_eq!(config.tile_size, DEFAULT_TILE_SIZE);
    assert_eq!(config.logging, LoggingArgs::default());
}

#[test]
fn ref_gc_bias_config_setters_update_programmatic_fields() {
    // Arrange
    let mut config = RefGCBiasConfig::new(
        PathBuf::from("ref.2bit"),
        PathBuf::from("out"),
        chromosomes(),
    );
    let updated_chromosomes = ChromosomeArgs {
        chromosomes: None,
        chromosomes_file: Some(PathBuf::from("chromosomes.txt")),
    };

    // Act
    config.set_ref_2bit(PathBuf::from("updated.2bit"));
    config.set_output_dir(PathBuf::from("updated_out"));
    config.set_output_prefix("hg38");
    config.set_n_threads(3);
    config.set_n_positions(12_345);
    config.set_seed(Some(99));
    config.set_windows(RefGCWindowsArgs {
        by_bed: Some(PathBuf::from("windows.bed")),
    });
    config.set_by_bed(Some(PathBuf::from("regions.bed")));
    config.set_chromosomes(updated_chromosomes.clone());
    config.set_blacklist(Some(vec![
        PathBuf::from("blacklist_a.bed"),
        PathBuf::from("blacklist_b.bed"),
    ]));
    config.set_fragment_lengths(FragmentLengthArgs {
        min_fragment_length: 20,
        max_fragment_length: 200,
    });
    config.fragment_lengths_mut().max_fragment_length = 220;
    config.set_end_offset(3);
    config.set_skip_interpolation(true);
    config.set_smoothing_sigma(1.25);
    config.set_smoothing_radius(4);
    config.set_skip_smoothing(true);
    config.set_tile_size(2_000_000);
    config.set_logging(LoggingArgs {
        log: LogSpec::Quiet,
    });

    // Assert
    assert_eq!(config.ref_genome.ref_2bit, PathBuf::from("updated.2bit"));
    assert_eq!(config.output_dir, PathBuf::from("updated_out"));
    assert_eq!(config.output_prefix, "hg38");
    assert_eq!(config.n_threads, 3);
    assert_eq!(config.n_positions, 12_345);
    assert_eq!(config.seed, Some(99));
    assert_eq!(
        config.windows,
        RefGCWindowsArgs {
            by_bed: Some(PathBuf::from("regions.bed")),
        }
    );
    assert_eq!(config.chromosomes, updated_chromosomes);
    assert_eq!(
        config.blacklist,
        Some(vec![
            PathBuf::from("blacklist_a.bed"),
            PathBuf::from("blacklist_b.bed"),
        ])
    );
    assert_eq!(
        config.fragment_lengths,
        FragmentLengthArgs {
            min_fragment_length: 20,
            max_fragment_length: 220,
        }
    );
    assert_eq!(config.end_offset, 3);
    assert!(config.skip_interpolation);
    assert_eq!(config.smoothing_sigma, 1.25);
    assert_eq!(config.smoothing_radius, 4);
    assert!(config.skip_smoothing);
    assert_eq!(config.tile_size, 2_000_000);
    assert_eq!(
        config.logging,
        LoggingArgs {
            log: LogSpec::Quiet,
        }
    );
}
