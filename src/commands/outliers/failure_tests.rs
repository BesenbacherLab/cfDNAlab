use super::*;
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs};

/// Construct failure-report settings without opening a BAM or changing global command state.
fn report_config(directory: &std::path::Path) -> OutliersConfig {
    OutliersConfig::new(
        IOCArgs {
            bam: "unopened_input.bam".into(),
            output_dir: directory.to_path_buf(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".into()]),
            chromosomes_file: None,
        },
    )
}

#[test]
fn contextual_failure_report_matches_companions_to_the_current_attempt() {
    let directory = tempfile::tempdir().unwrap();
    let mut config = report_config(directory.path());
    config.set_output_prefix("");
    let error = anyhow::Error::new(AnalysisFailures {
        stage: "model fitting",
        messages: vec!["chr1:0-500000 has no positive coverage".into()],
        histogram: Some("chromosome\tcoverage\tpositions\nchr1\t0\t500000\n".into()),
        audit: None,
    })
    .context("could not complete outlier detection");

    let report_path = write_failure_report(&config, &error).unwrap();

    assert_eq!(report_path, directory.path().join("outliers.error.txt"));
    let report = std::fs::read_to_string(report_path).unwrap();
    let histogram_path = directory.path().join("outliers.error_histogram.tsv");
    let histogram = std::fs::read_to_string(&histogram_path).unwrap();
    // Compare identifiers only within this attempt, without depending on wall-clock values
    let attempt = report
        .lines()
        .find_map(|line| line.strip_prefix("Attempt (Unix nanoseconds): "))
        .unwrap();
    assert!(histogram.starts_with(&format!(
        "# failed_attempt_unix_nanoseconds={attempt}\n# normal_results_withheld=true\n"
    )));
    assert!(histogram.ends_with("chromosome\tcoverage\tpositions\nchr1\t0\t500000\n"));
    assert!(report.contains("could not complete outlier detection"));
    assert!(report.contains("chr1:0-500000 has no positive coverage"));
    assert!(report.contains("unopened_input.bam"));
    assert!(report.contains(histogram_path.to_str().unwrap()));
    assert!(!directory.path().join("outliers.error_audit.tsv").exists());
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 2);
}

#[test]
fn retry_report_never_attributes_stale_companions_to_the_new_failure() {
    let directory = tempfile::tempdir().unwrap();
    let mut config = report_config(directory.path());
    config.set_output_prefix("sample");
    let old_histogram = directory.path().join("sample.outliers.error_histogram.tsv");
    let old_audit = directory.path().join("sample.outliers.error_audit.tsv");
    std::fs::write(&old_histogram, "previous attempt histogram").unwrap();
    std::fs::write(&old_audit, "previous attempt audit").unwrap();
    let analysis_error = anyhow::Error::new(AnalysisFailures {
        stage: "ambiguity policy",
        messages: vec!["chr1:10000-16000 exceeds the cutoff".into()],
        histogram: None,
        audit: Some("chromosome\tstart\tend\nchr1\t10000\t16000\n".into()),
    });

    let report_path = write_failure_report(&config, &analysis_error).unwrap();

    let report = std::fs::read_to_string(&report_path).unwrap();
    assert!(report.contains(old_audit.to_str().unwrap()));
    assert!(!report.contains(old_histogram.to_str().unwrap()));
    assert_eq!(
        std::fs::read_to_string(&old_histogram).unwrap(),
        "previous attempt histogram"
    );
    let audit = std::fs::read_to_string(&old_audit).unwrap();
    let attempt = report
        .lines()
        .find_map(|line| line.strip_prefix("Attempt (Unix nanoseconds): "))
        .unwrap();
    assert!(audit.starts_with(&format!("# failed_attempt_unix_nanoseconds={attempt}\n")));
    assert!(audit.ends_with("chr1\t10000\t16000\n"));

    // A subsequent failure without structured evidence must identify neither existing companion
    write_failure_report(&config, &anyhow::anyhow!("alignment reader failed")).unwrap();
    let report = std::fs::read_to_string(report_path).unwrap();
    assert!(report.contains("alignment reader failed"));
    assert!(report.contains("No context histogram or ambiguity audit was available"));
    assert!(!report.contains(old_audit.to_str().unwrap()));
    assert!(!report.contains(old_histogram.to_str().unwrap()));
    assert_eq!(std::fs::read_to_string(old_audit).unwrap(), audit);
}

#[test]
fn blocked_companion_reports_its_destination_without_touching_normal_results() {
    let directory = tempfile::tempdir().unwrap();
    let mut config = report_config(directory.path());
    config.set_output_prefix("sample");
    let models = directory.path().join("sample.outliers.models.tsv");
    let histogram = directory.path().join("sample.outliers.error_histogram.tsv");
    std::fs::write(&models, "previous normal results").unwrap();
    std::fs::create_dir(&histogram).unwrap();
    let error = anyhow::Error::new(AnalysisFailures {
        stage: "collected evidence",
        messages: vec!["current attempt failed".into()],
        histogram: Some("current histogram\n".into()),
        audit: Some("current audit\n".into()),
    });

    let failure = write_failure_report(&config, &error).unwrap_err();

    let diagnostic = format!("{failure:#}");
    assert!(diagnostic.contains("moving output file"), "{diagnostic}");
    assert!(
        diagnostic.contains(histogram.to_str().unwrap()),
        "{diagnostic}"
    );
    assert_eq!(
        std::fs::read_to_string(models).unwrap(),
        "previous normal results"
    );
    assert!(histogram.is_dir());
}
