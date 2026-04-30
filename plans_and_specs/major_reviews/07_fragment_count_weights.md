# `cfdna fragment-count-weights` review

Date: 2026-04-24

Scope: `src/commands/fragment_count_weights/*`, the shared smoothing-weight implementation in `src/commands/coverage_weights/*`, the internal `fcoverage` configuration boundary, downstream scaling-factor loading, README smoothing usage, and existing fragment-count-weight tests in `tests/test_normalize_genome_command.rs`, `tests/test_cli_smoke.rs`, and `src/commands/coverage_weights/coverage_weights_tests.rs`. I did not run tests.

Shared findings that affect this command:

- G-011 in `00_shared_package_notes.md`: scaling-factor TSV metadata is too thin for safe reuse.
- G-012 in `00_shared_package_notes.md`: short final stride bins are length-weighted only in the numerator.

## Release triage

Pre-release correctness/safety:

- G-012: short final stride bins distort edge scaling factors.
- G-011: scaling-factor TSVs need enough metadata for safe reuse.

Pre-release docs/schema polish:

- FCW-002: fragment-count outputs should not rely on coverage terminology.

## Findings

### FCW-002 - Low - Fragment-count outputs still use coverage terminology

The command documentation correctly frames the raw values as fragment-count density/support rather than ordinary coverage ([config.rs](../../src/commands/fragment_count_weights/config.rs#L32-L39), [config.rs](../../src/commands/fragment_count_weights/config.rs#L57-L60)). The shared writer still logs "overlapping position-coverage" and writes columns named `average_pos_coverage` and `average_overlapping_pos_coverage` for the fragment-count output ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L139-L142), [coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L159-L163)). The CLI smoke test also pins those coverage-named columns for `fragment-count-weights` ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L366-L375)).

Impact: the file name says `fragment_counts`, but the primary numeric columns read like coverage-weight output. Since downstream loaders use only `scaling_factor`, this is mostly a human-facing schema clarity issue, but it makes it easier to miss whether a TSV was built from coverage support or fragment-count support.

Recommended fix:

- Rename the shared columns to neutral names like `average_pos_support` and `average_overlapping_pos_support`, or emit command-specific column names while keeping `scaling_factor` stable.
- Add source metadata as part of G-011 so downstream tools do not need to infer the smoothing source from filenames or column wording.
- Update the CLI smoke and parsing tests to pin the intended fragment-count schema.

## Existing Coverage Notes

The command already has a basic integration check that it writes one stride row per chromosome span ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L260-L290)), an identity-smoothing comparison proving fragment-count mode differs from coverage mode for mixed fragment lengths ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L343-L480)), a CLI smoke test for the output path, metadata line, header, and row count ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L319-L383)), and a helper test proving the shared builder uses `LengthNormalizationMode::UnitMass` for fragment-count weights ([coverage_weights_tests.rs](../../src/commands/coverage_weights/coverage_weights_tests.rs#L395-L416)).

The important missing coverage from this review is fragment-count expected values with real smoothing (`bin_size > stride`), short-final-bin edge behavior from G-012, richer scaling metadata from G-011, and the intended fragment-count output terminology.
