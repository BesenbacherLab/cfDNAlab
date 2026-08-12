use super::{
    ChromosomeOutlierWeights, load_outlier_weights_tsv, minimum_overlapping_keep_weight,
};
use crate::shared::bam::Contigs;
use crate::shared::interval::{IndexedInterval, Interval};
use crate::shared::overlaps::find_overlapping_windows;
use fxhash::FxHashMap;
use std::io::Write;
use tempfile::NamedTempFile;

fn selected_contigs() -> Contigs {
    let mut contigs = FxHashMap::with_hasher(Default::default());
    contigs.insert("chr1".to_string(), (0, 100));
    contigs.insert("chr2".to_string(), (1, 80));
    Contigs { contigs }
}

fn write_input(contents: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("temporary outlier input should be created");
    file.write_all(contents.as_bytes())
        .expect("temporary outlier input should be written");
    file
}

#[test]
fn loader_accepts_sparse_output_metadata_and_tracks_the_minimum_positive_weight() {
    // Arrange
    // chr2 has no rows, which means its fragments retain the implicit weight 1.0. The zero row
    // can exclude fragments but cannot define the smallest positive support used by cleanup.
    let file = write_input(
        "# omitted_keep_weight=1.0\n\
         # fragment_overlap_rule=minimum_keep_weight\n\
         CHROMOSOME\tSTART\tEND\tKEEP_WEIGHT\n\
         chr1\t10\t20\t0.5\n\
         chr1\t30\t40\t0\n",
    );
    let chromosomes = vec!["chr1".to_string(), "chr2".to_string()];

    // Act
    let loaded = load_outlier_weights_tsv(file.path(), &chromosomes, &selected_contigs())
        .expect("valid sparse outlier weights should load");

    // Assert
    let chr1 = loaded
        .by_chromosome
        .get("chr1")
        .expect("chr1 weights should be present");
    assert_eq!(chr1.intervals.len(), 2);
    assert_eq!(chr1.keep_weights, vec![0.5, 0.0]);
    assert!(!loaded.by_chromosome.contains_key("chr2"));
    assert_eq!(loaded.minimum_positive_keep_weight, 0.5);
}

#[test]
fn zero_keep_weights_do_not_lower_the_minimum_positive_weight() {
    // Arrange
    // A zero weight removes overlapping fragment contributions. It is not positive support and
    // therefore must not lower the cleanup bound used for remaining implicit weight 1.0 regions.
    let file = write_input(
        "chromosome\tstart\tend\tkeep_weight\n\
         chr1\t10\t20\t0\n",
    );

    // Act
    let loaded = load_outlier_weights_tsv(file.path(), &["chr1".to_string()], &selected_contigs())
        .expect("zero keep weights are valid");

    // Assert
    assert_eq!(loaded.minimum_positive_keep_weight, 1.0);
}

#[test]
fn overlap_resolution_uses_any_positive_overlap_and_the_smallest_regional_weight() {
    // Arrange
    // The fragment [19, 31) overlaps each outlier interval by exactly 1 bp. Full regional weights
    // apply despite the small overlaps, and the smaller of the two weights should win.
    let chromosome_weights = ChromosomeOutlierWeights {
        intervals: vec![
            IndexedInterval::new(10, 20, 0).expect("valid interval"),
            IndexedInterval::new(30, 40, 1).expect("valid interval"),
        ],
        keep_weights: vec![0.6, 0.2],
    };
    let fragment = Interval::new(19, 31).expect("valid fragment interval");
    let mut interval_pointer = 0usize;

    // Act
    let overlaps = find_overlapping_windows(
        100,
        &mut interval_pointer,
        Some(&chromosome_weights.intervals),
        None,
        fragment,
        1.0 / 101.0,
        100,
    )
    .expect("overlap sweep should succeed");
    let keep_weight = minimum_overlapping_keep_weight(overlaps.as_ref(), &chromosome_weights)
        .expect("overlap weights should resolve");

    // Assert
    assert_eq!(keep_weight, 0.2);
}

#[test]
fn overlap_resolution_treats_touching_half_open_intervals_as_no_overlap() {
    // Arrange
    let chromosome_weights = ChromosomeOutlierWeights {
        intervals: vec![IndexedInterval::new(10, 20, 0).expect("valid interval")],
        keep_weights: vec![0.1],
    };
    let touching_fragment = Interval::new(20, 30).expect("valid fragment interval");
    let mut interval_pointer = 0usize;

    // Act
    let overlaps = find_overlapping_windows(
        100,
        &mut interval_pointer,
        Some(&chromosome_weights.intervals),
        None,
        touching_fragment,
        1.0 / 101.0,
        100,
    )
    .expect("overlap sweep should succeed");
    let keep_weight = minimum_overlapping_keep_weight(overlaps.as_ref(), &chromosome_weights)
        .expect("missing overlap should resolve to identity");

    // Assert
    assert!(overlaps.is_none());
    assert_eq!(keep_weight, 1.0);
}

#[test]
fn loader_rejects_unsorted_or_overlapping_rows() {
    // Arrange
    let file = write_input(
        "chromosome\tstart\tend\tkeep_weight\n\
         chr1\t30\t40\t0.5\n\
         chr1\t20\t35\t0.2\n",
    );

    // Act
    let error = load_outlier_weights_tsv(file.path(), &["chr1".to_string()], &selected_contigs())
        .expect_err("unsorted overlapping rows should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("must be sorted and non-overlapping"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn loader_rejects_non_finite_or_out_of_range_keep_weights() {
    for invalid_weight in ["NaN", "-0.1", "1.1"] {
        // Arrange
        let file = write_input(&format!(
            "chromosome\tstart\tend\tkeep_weight\nchr1\t10\t20\t{invalid_weight}\n"
        ));

        // Act
        let error =
            load_outlier_weights_tsv(file.path(), &["chr1".to_string()], &selected_contigs())
                .expect_err("invalid keep weight should fail");

        // Assert
        assert!(
            error
                .to_string()
                .contains("keep_weight must be finite and between 0.0 and 1.0 inclusive"),
            "unexpected error for {invalid_weight}: {error:#}"
        );
    }
}

#[test]
fn loader_rejects_intervals_beyond_the_selected_bam_contig() {
    // Arrange
    let file = write_input(
        "chromosome\tstart\tend\tkeep_weight\n\
         chr1\t90\t101\t0.5\n",
    );

    // Act
    let error = load_outlier_weights_tsv(file.path(), &["chr1".to_string()], &selected_contigs())
        .expect_err("out-of-contig interval should fail");

    // Assert
    assert!(
        error.to_string().contains("beyond chromosome length 100"),
        "unexpected error: {error:#}"
    );
}
