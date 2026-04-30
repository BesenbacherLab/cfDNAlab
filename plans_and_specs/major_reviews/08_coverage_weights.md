# `cfdna coverage-weights` review

Date: 2026-04-24

Scope: `src/commands/coverage_weights/*`, the internal `fcoverage` configuration boundary, downstream scaling-factor loading in `src/shared/scale_genome.rs`, README smoothing usage, and existing coverage-weight tests in `tests/test_normalize_genome_command.rs`, `tests/test_cli_smoke.rs`, `tests/test_cross_command_artifact_matrix.rs`, and `src/commands/coverage_weights/coverage_weights_tests.rs`. I did not run tests.

Shared findings that affect this command:

- G-011 in `00_shared_package_notes.md`: scaling-factor TSV metadata is too thin for safe reuse.
- G-012 in `00_shared_package_notes.md`: short final stride bins are length-weighted only in the numerator.

## Release triage

Pre-release correctness/safety:

- G-012: short final stride bins distort edge scaling factors.
- G-011: scaling-factor TSVs need enough metadata for safe reuse.

## Findings

No additional `coverage-weights`-specific findings beyond the shared smoothing-weight findings above.

The command's feature wiring is internally consistent in the current feature matrix: `cmd_coverage_weights` enables `cmd_fcoverage` ([Cargo.toml](../../Cargo.toml#L36-L39)), the module is gated by `cmd_coverage_weights` ([mod.rs](../../src/commands/mod.rs#L7-L8)), and the CLI subcommand is gated by the same feature ([cli_app.rs](../../src/cli_app.rs#L5-L6), [cli_app.rs](../../src/cli_app.rs#L79-L80)). The command is also a thin wrapper over the shared smoothing implementation, with coverage mode selected by `normalize_by_length = false` ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L89-L90), [coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L233-L238)).

## Existing Coverage Notes

The command already has direct coverage for output row ranges and non-zero support behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L208-L258)), contiguous output bins ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L293-L341)), hand-derived smoothing and scaling values for a simple fragment ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L484-L563)), unpaired read-as-fragment behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L565-L640)), bin-size validation ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L642-L689)), support-floor and length-weighted global mean behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L691-L907)), one shared global mean across chromosomes ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L909-L1059)), default MAPQ behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L1061-L1274)), chromosome row ordering ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L1276-L1373)), blacklist masking ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L1375-L1457)), CLI smoke output shape ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L253-L317)), and cross-command consumption of a real scaling TSV with GC correction ([test_cross_command_artifact_matrix.rs](../../tests/test_cross_command_artifact_matrix.rs#L96-L147), [test_cross_command_artifact_matrix.rs](../../tests/test_cross_command_artifact_matrix.rs#L240-L450)).

The important missing coverage from this review is standalone feature-matrix compilation, short-final-bin smoothing from G-012, and richer scaling metadata from G-011.
