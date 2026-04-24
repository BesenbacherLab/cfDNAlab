#![cfg(feature = "cli")]

mod fixtures;

use anyhow::{Context, Result, bail};
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

fn binary_name(base_name: &str) -> String {
    if cfg!(windows) {
        format!("{base_name}.exe")
    } else {
        base_name.to_string()
    }
}

fn cfdna_bin_path() -> Result<PathBuf> {
    // Preferred path when Cargo exports the integration-test binary location.
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_cfdna") {
        return Ok(PathBuf::from(path));
    }

    // Fallback for IDE runners that do not export CARGO_BIN_EXE_*.
    // We try target/{profile}/cfdna first, then target/{profile}/deps/cfdna-<hash>.
    let current_exe = std::env::current_exe().context("failed to read current test binary path")?;
    let deps_dir = current_exe
        .parent()
        .context("failed to derive deps directory from current test binary path")?;
    let target_dir = deps_dir
        .parent()
        .context("failed to derive target directory from deps path")?;

    let direct_path = target_dir.join(binary_name("cfdna"));
    if direct_path.is_file() {
        return Ok(direct_path);
    }

    let mut hashed_candidates = Vec::new();
    for entry in std::fs::read_dir(deps_dir).with_context(|| {
        format!(
            "failed to list candidate binaries in {}",
            deps_dir.display()
        )
    })? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let file_name = match path.file_name().and_then(OsStr::to_str) {
            Some(name) => name,
            None => continue,
        };
        let extension = path.extension().and_then(OsStr::to_str);
        let looks_like_hashed_binary = file_name.starts_with("cfdna-");
        let is_makefile_dep = extension == Some("d");
        if looks_like_hashed_binary && !is_makefile_dep {
            hashed_candidates.push(path);
        }
    }
    hashed_candidates.sort_by_key(|path| {
        std::fs::metadata(path)
            .and_then(|meta| meta.modified())
            .ok()
            .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
    });
    if let Some(path) = hashed_candidates.into_iter().last() {
        return Ok(path);
    }

    bail!(
        "Could not locate cfdna binary. Tried CARGO_BIN_EXE_cfdna, {}, and hashed binaries under {}",
        direct_path.display(),
        deps_dir.display()
    );
}

fn command_output(command_name: &str, args: &[&str]) -> Result<std::process::Output> {
    let output = Command::new(cfdna_bin_path()?)
        .arg(command_name)
        .args(args)
        .output()
        .with_context(|| format!("failed running cfdna {command_name} {}", args.join(" ")))?;
    Ok(output)
}

fn assert_success_with_logs(output: &std::process::Output, command_desc: &str) {
    let stdout_text = String::from_utf8_lossy(&output.stdout);
    let stderr_text = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "Expected {command_desc} to succeed.\nstdout:\n{stdout_text}\nstderr:\n{stderr_text}"
    );
}

fn path_text(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

#[test]
fn help_text_is_available_for_all_enabled_release_commands() -> Result<()> {
    // Human verification status: unverified
    let mut release_commands = Vec::new();
    #[cfg(feature = "cmd_gc_bias")]
    release_commands.push("gc-bias");
    #[cfg(feature = "cmd_ref_gc_bias")]
    release_commands.push("ref-gc-bias");
    #[cfg(feature = "cmd_coverage_weights")]
    release_commands.push("coverage-weights");
    #[cfg(feature = "cmd_fragment_count_weights")]
    release_commands.push("fragment-count-weights");
    #[cfg(feature = "cmd_lengths")]
    release_commands.push("lengths");
    #[cfg(feature = "cmd_fcoverage")]
    release_commands.push("fcoverage");
    #[cfg(feature = "cmd_midpoints")]
    release_commands.push("midpoints");
    #[cfg(feature = "cmd_bam_to_bam")]
    release_commands.push("bam-to-bam");
    #[cfg(feature = "cmd_bam_to_frag")]
    release_commands.push("bam-to-frag");
    #[cfg(feature = "cmd_frag_to_bam")]
    release_commands.push("frag-to-bam");

    assert!(
        !release_commands.is_empty(),
        "Expected at least one release command to be enabled in this build"
    );

    for command_name in release_commands {
        let output = command_output(command_name, &["--help"])?;
        let stdout_text = String::from_utf8_lossy(&output.stdout);
        let stderr_text = String::from_utf8_lossy(&output.stderr);

        assert!(
            output.status.success(),
            "Expected cfdna {command_name} --help to succeed.\nstdout:\n{stdout_text}\nstderr:\n{stderr_text}"
        );
        assert!(
            !stdout_text.trim().is_empty(),
            "Expected cfdna {command_name} --help to print text output"
        );
        assert!(
            stdout_text.contains("Usage:"),
            "Expected help output for cfdna {command_name} to contain a Usage section"
        );
    }

    Ok(())
}

#[cfg(feature = "cmd_ends")]
#[test]
fn ends_help_only_shows_collapse_complements_when_experimental_feature_is_enabled() -> Result<()> {
    // Human verification status: unverified
    let output = command_output("ends", &["--help"])?;
    let stdout_text = String::from_utf8_lossy(&output.stdout);
    let stderr_text = String::from_utf8_lossy(&output.stderr);

    assert!(
        output.status.success(),
        "Expected cfdna ends --help to succeed.\nstdout:\n{stdout_text}\nstderr:\n{stderr_text}"
    );

    if cfg!(feature = "ends_experimental") {
        assert!(
            stdout_text.contains("--collapse-complement"),
            "Expected experimental ends help to show --collapse-complement.\nstdout:\n{stdout_text}"
        );
    } else {
        assert!(
            !stdout_text.contains("--collapse-complement"),
            "Expected default ends help to hide --collapse-complement.\nstdout:\n{stdout_text}"
        );
    }

    Ok(())
}

#[cfg(feature = "cmd_lengths")]
#[test]
fn lengths_cli_minimal_invocation_writes_output_files_with_expected_prefix() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // The command contract says lengths writes:
    // - <prefix>.length_counts.npy
    // - <prefix>.fragment_length_settings.json
    // We run a minimal binary invocation with a tiny deterministic BAM fixture.
    let bam_fixture = fixtures::simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let output_prefix = "cli_smoke_lengths";

    // Act
    let bam_path = path_text(&bam_fixture.bam);
    let output_path = path_text(out_dir.path());
    let output = command_output(
        "lengths",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            output_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            output_prefix,
        ],
    )?;

    assert_success_with_logs(&output, "cfdna lengths minimal invocation");

    // Assert
    let npy_path = out_dir
        .path()
        .join(format!("{output_prefix}.length_counts.npy"));
    let settings_path = out_dir
        .path()
        .join(format!("{output_prefix}.fragment_length_settings.json"));
    assert!(
        npy_path.exists(),
        "Expected output file to exist: {}",
        npy_path.display()
    );
    assert!(
        settings_path.exists(),
        "Expected output file to exist: {}",
        settings_path.display()
    );

    let settings_text = std::fs::read_to_string(&settings_path)?;
    assert!(
        settings_text.contains("\"min_fragment_length\":10"),
        "Expected settings file to include min fragment length"
    );
    assert!(
        settings_text.contains("\"max_fragment_length\":200"),
        "Expected settings file to include max fragment length"
    );

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn coverage_weights_cli_minimal_invocation_writes_scaling_tsv() -> Result<()> {
    // Human verification status: unverified
    // Arrange: simple_inward_bam has chr1 length 200 and one fragment spanning [20,80).
    // With stride 20 this yields exactly 10 stride bins -> 13 TSV lines including two metadata
    // lines and one header line.
    let bam_fixture = fixtures::simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let output_prefix = "cli_smoke_covweights";
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "coverage-weights",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--bin-size",
            "40",
            "--stride",
            "20",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            output_prefix,
        ],
    )?;
    assert_success_with_logs(&output, "cfdna coverage-weights minimal invocation");

    // Assert
    let scaling_path = out_dir
        .path()
        .join(format!("{output_prefix}.coverage.scaling_factors.tsv"));
    assert!(scaling_path.exists(), "Expected {}", scaling_path.display());

    let content = std::fs::read_to_string(&scaling_path)?;
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(
        lines.first().copied().unwrap_or_default(),
        "# gc_mode=uncorrected"
    );
    assert_eq!(
        lines.get(1).copied().unwrap_or_default(),
        "# ignore_gap=false"
    );
    assert_eq!(
        lines.get(2).copied().unwrap_or_default(),
        "chromosome\tstart\tend\taverage_pos_coverage\taverage_overlapping_pos_coverage\tscaling_factor"
    );
    assert_eq!(
        lines.len(),
        13,
        "Expected 2 metadata lines + 1 header + 10 stride-bin rows for chr1 length 200 with stride 20"
    );

    Ok(())
}

#[cfg(feature = "cmd_fragment_count_weights")]
#[test]
fn fragment_count_weights_cli_minimal_invocation_writes_scaling_tsv() -> Result<()> {
    // Human verification status: unverified
    // Arrange: simple_inward_bam has chr1 length 200 and one fragment spanning [20,80).
    // With stride 20 this yields exactly 10 stride bins -> 12 TSV lines including one metadata
    // line and one header line.
    let bam_fixture = fixtures::simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let output_prefix = "cli_smoke_fragcountweights";
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "fragment-count-weights",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--bin-size",
            "40",
            "--stride",
            "20",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            output_prefix,
        ],
    )?;
    assert_success_with_logs(&output, "cfdna fragment-count-weights minimal invocation");

    // Assert
    let scaling_path = out_dir.path().join(format!(
        "{output_prefix}.fragment_counts.scaling_factors.tsv"
    ));
    assert!(scaling_path.exists(), "Expected {}", scaling_path.display());

    let content = std::fs::read_to_string(&scaling_path)?;
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(
        lines.first().copied().unwrap_or_default(),
        "# gc_mode=uncorrected"
    );
    assert_eq!(
        lines.get(1).copied().unwrap_or_default(),
        "chromosome\tstart\tend\taverage_pos_coverage\taverage_overlapping_pos_coverage\tscaling_factor"
    );
    assert_eq!(
        lines.len(),
        12,
        "Expected 1 metadata line + 1 header + 10 stride-bin rows for chr1 length 200 with stride 20"
    );

    Ok(())
}

#[cfg(feature = "cmd_fcoverage")]
#[test]
fn fcoverage_cli_minimal_invocation_writes_expected_positional_run() -> Result<()> {
    // Human verification status: unverified
    // Arrange: simple_inward_bam has one fragment spanning [20,80) on chr1.
    // In plain positional mode without correction, expected run is coverage 1 on [20,80).
    let bam_fixture = fixtures::simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let output_prefix = "cli_smoke_fcoverage";
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "fcoverage",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            output_prefix,
        ],
    )?;
    assert_success_with_logs(&output, "cfdna fcoverage minimal invocation");

    // Assert
    let coverage_path = out_dir.path().join(format!(
        "{output_prefix}.fcoverage.per_position.bedgraph.zst"
    ));
    assert!(
        coverage_path.exists(),
        "Expected {}",
        coverage_path.display()
    );

    let coverage_text = fixtures::read_zst_to_string(&coverage_path)?;
    assert!(
        coverage_text.contains("chr1\t20\t80\t1"),
        "Expected a single-fragment positional run in output"
    );

    Ok(())
}

#[cfg(feature = "cmd_midpoints")]
#[test]
fn midpoints_cli_minimal_invocation_writes_profiles_and_group_index() -> Result<()> {
    // Human verification status: unverified
    // Arrange: one window in one group with one fragment-length bin.
    let bam_fixture = fixtures::simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let intervals_path = out_dir.path().join("intervals.bed");
    fixtures::write_bed(&intervals_path, &[("chr1", 20, 60, "groupA")])?;
    let output_prefix = "cli_smoke_midpoints";
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());
    let intervals_text = path_text(&intervals_path);

    // Act
    let output = command_output(
        "midpoints",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--intervals",
            intervals_text.as_str(),
            "--min-mapq",
            "0",
            "--length-bins",
            "20",
            "120",
            "--output-prefix",
            output_prefix,
        ],
    )?;
    assert_success_with_logs(&output, "cfdna midpoints minimal invocation");

    // Assert
    let profiles_path = out_dir
        .path()
        .join(format!("{output_prefix}.midpoint_profiles.npy"));
    let group_index_path = out_dir
        .path()
        .join(format!("{output_prefix}.group_index.tsv"));
    assert!(
        profiles_path.exists(),
        "Expected {}",
        profiles_path.display()
    );
    assert!(
        group_index_path.exists(),
        "Expected {}",
        group_index_path.display()
    );
    let group_index_text = std::fs::read_to_string(&group_index_path)?;
    assert!(group_index_text.contains("groupA"));

    Ok(())
}

#[cfg(feature = "cmd_bam_to_bam")]
#[test]
fn bam_to_bam_cli_minimal_invocation_writes_output_bam() -> Result<()> {
    // Human verification status: unverified
    // Arrange
    let bam_fixture = fixtures::simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let out_bam = out_dir.path().join("smoke.bam");
    let bam_path = path_text(&bam_fixture.bam);
    let out_bam_text = path_text(&out_bam);

    // Act
    let output = command_output(
        "bam-to-bam",
        &[
            "--in-bam",
            bam_path.as_str(),
            "--out-bam",
            out_bam_text.as_str(),
            "--chromosomes",
            "chr1",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
        ],
    )?;
    assert_success_with_logs(&output, "cfdna bam-to-bam minimal invocation");

    // Assert
    assert!(out_bam.exists(), "Expected {}", out_bam.display());
    let file_size = std::fs::metadata(&out_bam)?.len();
    assert!(file_size > 0, "Expected non-empty BAM output");

    Ok(())
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn bam_to_frag_cli_minimal_invocation_writes_frag_and_header_files() -> Result<()> {
    // Human verification status: unverified
    // Arrange
    let bam_fixture = fixtures::simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let output_prefix = "cli_smoke_fragments";
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "bam-to-frag",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            output_prefix,
        ],
    )?;
    assert_success_with_logs(&output, "cfdna bam-to-frag minimal invocation");

    // Assert
    let frag_path = out_dir.path().join(format!("{output_prefix}.frag.tsv.gz"));
    let header_path = out_dir
        .path()
        .join(format!("{output_prefix}.frag.header.tsv"));
    assert!(frag_path.exists(), "Expected {}", frag_path.display());
    assert!(header_path.exists(), "Expected {}", header_path.display());

    Ok(())
}

#[cfg(feature = "cmd_frag_to_bam")]
#[test]
fn frag_to_bam_cli_minimal_invocation_writes_output_bam() -> Result<()> {
    // Human verification status: unverified
    // Arrange: one valid frag row and one matching chrom.sizes entry.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    std::fs::write(&frag_path, "chr1\t20\t80\t60\t+\n")?;
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");
    std::fs::write(&chrom_sizes_path, "chr1\t200\n")?;
    let output_prefix = "cli_smoke_frag_to_bam";
    let frag_text = path_text(&frag_path);
    let chrom_sizes_text = path_text(&chrom_sizes_path);
    let output_dir_text = path_text(output_dir.path());

    // Act
    let output = command_output(
        "frag-to-bam",
        &[
            "--frag",
            frag_text.as_str(),
            "--output-dir",
            output_dir_text.as_str(),
            "--chrom-sizes",
            chrom_sizes_text.as_str(),
            "--chromosomes",
            "chr1",
            "--output-prefix",
            output_prefix,
            "--min-mapq",
            "0",
        ],
    )?;
    assert_success_with_logs(&output, "cfdna frag-to-bam minimal invocation");

    // Assert
    let bam_path = output_dir
        .path()
        .join(format!("{output_prefix}.fragments.bam"));
    assert!(bam_path.exists(), "Expected {}", bam_path.display());
    let file_size = std::fs::metadata(&bam_path)?.len();
    assert!(file_size > 0, "Expected non-empty BAM output");

    Ok(())
}

#[cfg(feature = "cmd_ref_gc_bias")]
#[test]
fn ref_gc_bias_cli_minimal_invocation_writes_reference_package() -> Result<()> {
    // Human verification status: unverified
    // Arrange: Use tiny deterministic reference and conservative settings.
    // With `--output-prefix`, the command contract says the package should be written as
    // `<prefix>.ref_gc_package.npz`.
    let reference = fixtures::simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let output_prefix = "cli_smoke_ref_gc";
    let ref_path = path_text(&reference.path);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "ref-gc-bias",
        &[
            "--ref-2bit",
            ref_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--output-prefix",
            output_prefix,
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--n-positions",
            "100",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "120",
            "--end-offset",
            "0",
            "--tile-size",
            "1000000",
        ],
    )?;
    assert_success_with_logs(&output, "cfdna ref-gc-bias minimal invocation");

    // Assert
    let package_path = out_dir
        .path()
        .join(format!("{output_prefix}.ref_gc_package.npz"));
    assert!(package_path.exists(), "Expected {}", package_path.display());
    assert!(
        !out_dir.path().join("ref_gc_package.npz").exists(),
        "Did not expect unprefixed package when --output-prefix is supplied"
    );

    Ok(())
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
#[test]
fn gc_bias_cli_minimal_invocation_writes_correction_package() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // 1) Build reference package from tiny deterministic reference.
    // 2) Run gc-bias on tiny deterministic BAM against that reference package.
    let bam_fixture = fixtures::simple_inward_bam()?;
    let reference = fixtures::simple_reference_twobit()?;
    let ref_gc_dir = TempDir::new()?;
    let gc_out_dir = TempDir::new()?;
    let ref_gc_prefix = "cli_smoke_ref_gc";

    let ref_path = path_text(&reference.path);
    let ref_gc_out = path_text(ref_gc_dir.path());
    let ref_gc_file = path_text(
        &ref_gc_dir
            .path()
            .join(format!("{ref_gc_prefix}.ref_gc_package.npz")),
    );
    let ref_gc_output = command_output(
        "ref-gc-bias",
        &[
            "--ref-2bit",
            ref_path.as_str(),
            "--output-dir",
            ref_gc_out.as_str(),
            "--output-prefix",
            ref_gc_prefix,
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--n-positions",
            "100",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "120",
            "--end-offset",
            "0",
            "--tile-size",
            "1000000",
        ],
    )?;
    assert_success_with_logs(
        &ref_gc_output,
        "cfdna ref-gc-bias pre-step for gc-bias smoke test",
    );
    assert!(
        ref_gc_dir
            .path()
            .join(format!("{ref_gc_prefix}.ref_gc_package.npz"))
            .exists(),
        "Expected reference package before running gc-bias"
    );

    // Act
    let bam_path = path_text(&bam_fixture.bam);
    let gc_out_path = path_text(gc_out_dir.path());
    let gc_output = command_output(
        "gc-bias",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            gc_out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--ref-2bit",
            ref_path.as_str(),
            "--ref-gc-file",
            ref_gc_file.as_str(),
            "--global",
            "--min-mapq",
            "0",
        ],
    )?;
    assert_success_with_logs(&gc_output, "cfdna gc-bias minimal invocation");

    // Assert
    let correction_path = gc_out_dir.path().join("gc_bias_correction.npz");
    assert!(
        correction_path.exists(),
        "Expected {}",
        correction_path.display()
    );

    Ok(())
}
