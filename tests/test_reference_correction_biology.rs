#![cfg(all(feature = "cmd_ends", feature = "cmd_ref_kmers", feature = "testing"))]

use anyhow::Result;
use cfdnalab::{
    RunOptions,
    output_loaders::{load_ends_output, load_ref_kmers_output},
    run_like_cli::{
        common::{ChromosomeArgs, IOCArgs},
        ends::{ClipStrategy, EndsConfig, KmerSource, run_ends},
        ref_kmers::{RefKmerOrientation, RefKmersConfig, run_ref_kmers},
    },
    testing::{PairedFragmentSpec, bam_from_fragments, twobit_from_sequences},
};
use tempfile::TempDir;

fn assert_close(observed: f64, expected: f64) {
    assert!(
        (observed - expected).abs() < 1e-12,
        "observed {observed}, expected {expected}"
    );
}

/// Verify reference correction removes the composition bias from fragments sampled in proportion
/// to the available reference end motifs.
#[test]
fn unbiased_reference_end_sampling_has_uniform_corrected_motif_counts() -> Result<()> {
    // Arrange:
    // - The 81 bp reference is AACG repeated 20 times followed by A. Its 80 dinucleotide starts
    //   contain AA=20, AC=20, CG=20, and GA=20.
    // - `orientation=both` shares non-palindromic counts with their reverse complements. The
    //   resulting positive counts are AA=10, AC=10, CG=20, GA=10, GT=10, TC=10, and TT=10.
    // - A 14 bp fragment has its left and raw right dinucleotides 12 bp apart. Since 12 is
    //   divisible by the 4 bp repeat length, they are the same. One fragment at every start from
    //   0 through 67 samples each repeat phase 17 times. Right-end orientation therefore gives
    //   sample counts AA=17, AC=17, CG=34, GA=17, GT=17, TC=17, and TT=17.
    // - The 136 observed ends exactly follow the reference frequencies. Seven supported motifs
    //   therefore all correct to 136/7 while the total remains 136. The non-palindromic pairs
    //   AC/GT and GA/TC also make this sensitive to reversal order at the right end.
    let reference = twobit_from_sequences(
        "reference_correction_biology",
        vec![("chr1".to_string(), format!("{}A", "AACG".repeat(20)))],
    )?;
    let fragments = (0..68)
        .map(|fragment_start| PairedFragmentSpec::new(0, fragment_start, 14, 4).build())
        .collect::<Result<Vec<_>>>()?;
    let bam = bam_from_fragments(
        "reference_correction_biology",
        vec![("chr1".to_string(), 81)],
        fragments,
        Vec::new(),
    )?;
    let output_dir = TempDir::new()?;
    let chromosomes = ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string()]),
        chromosomes_file: None,
    };

    let mut ends_config = EndsConfig::new(
        IOCArgs {
            bam: bam.bam,
            output_dir: output_dir.path().to_path_buf(),
            n_threads: 1,
        },
        chromosomes.clone(),
        2,
        0,
    );
    ends_config.output_prefix = "unbiased".to_string();
    ends_config.set_ref_2bit(Some(reference.path.clone()));
    ends_config.source_inside = KmerSource::Reference;
    ends_config.all_motifs = true;
    ends_config.clip.clip_strategy = ClipStrategy::Aligned;
    ends_config.set_min_mapq(0);
    ends_config.set_require_proper_pair(false);
    ends_config.set_tile_size(1_000_000);
    {
        let fragment_lengths = ends_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 14;
        fragment_lengths.max_fragment_length = 14;
    }
    run_ends(&ends_config, RunOptions::new_quiet())?;

    let mut ref_kmers_config = RefKmersConfig::new(
        reference.path,
        output_dir.path().to_path_buf(),
        2,
        chromosomes,
    );
    ref_kmers_config.set_output_prefix("unbiased");
    ref_kmers_config.set_n_threads(1);
    ref_kmers_config.set_orientation(RefKmerOrientation::Both);
    ref_kmers_config.set_all_motifs(true);
    ref_kmers_config.set_tile_size(1_000_000);
    run_ref_kmers(&ref_kmers_config, RunOptions::new_quiet())?;

    // Act
    let ends = load_ends_output(output_dir.path().join("unbiased.end_motifs.zarr"))?;
    let ref_kmers = load_ref_kmers_output(output_dir.path().join("unbiased.ref_kmers.zarr"))?;
    for (motif, expected_count, expected_reference_frequency) in [
        ("AA", 17.0, 1.0 / 8.0),
        ("AC", 17.0, 1.0 / 8.0),
        ("CG", 34.0, 1.0 / 4.0),
        ("GA", 17.0, 1.0 / 8.0),
        ("GT", 17.0, 1.0 / 8.0),
        ("TC", 17.0, 1.0 / 8.0),
        ("TT", 17.0, 1.0 / 8.0),
    ] {
        assert_close(
            ends.count_for_motif(0, &format!("_{motif}"))?.unwrap(),
            expected_count,
        );
        assert_close(
            ref_kmers.frequency_for_motif(0, motif)?.unwrap(),
            expected_reference_frequency,
        );
    }
    let corrected = ends.select_corrected_counts(&ref_kmers).read()?;

    // Assert
    assert_eq!(corrected.shape(), (1, 16));
    for motif in ["_AA", "_AC", "_CG", "_GA", "_GT", "_TC", "_TT"] {
        let motif_index = corrected
            .motif_labels()
            .iter()
            .position(|label| label == motif)
            .expect("supported motif should be present");
        assert_close(corrected.count(0, motif_index).unwrap(), 136.0 / 7.0);
    }
    assert_close(
        corrected.to_dense_matrix()?.values_row_major().iter().sum(),
        136.0,
    );
    Ok(())
}
