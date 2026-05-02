# `cfdna fragment-count-weights` review

Date: 2026-04-24

Scope: `src/commands/fragment_count_weights/*`, the shared smoothing-weight implementation in `src/commands/coverage_weights/*`, the internal `fcoverage` configuration boundary, downstream scaling-factor loading, README smoothing usage, and existing fragment-count-weight tests in `tests/test_normalize_genome_command.rs`, `tests/test_cli_smoke.rs`, and `src/commands/coverage_weights/coverage_weights_tests.rs`. I did not run tests.

Shared findings that affect this command:

- No active shared correctness findings from `00_shared_package_notes.md`; remaining items below are command-specific.

## Release triage

Pre-release correctness/safety:

- None currently active.

## Findings

No active findings.

## Existing Coverage Notes

The command already has a basic integration check that it writes one stride row per chromosome span ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L260-L290)), an identity-smoothing comparison proving fragment-count mode differs from coverage mode for mixed fragment lengths ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L343-L480)), a CLI smoke test for the output path, metadata line, header, and row count ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L319-L383)), and a helper test proving the shared builder uses `LengthNormalizationMode::UnitMass` for fragment-count weights ([coverage_weights_tests.rs](../../src/commands/coverage_weights/coverage_weights_tests.rs#L395-L416)).

The important missing coverage from this review is fragment-count expected values with real smoothing (`bin_size > stride`).
