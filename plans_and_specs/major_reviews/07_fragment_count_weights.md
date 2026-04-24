# `cfdna fragment-count-weights` review

Date: 2026-04-24

Scope: `src/commands/fragment_count_weights/*`, the shared smoothing-weight implementation in `src/commands/coverage_weights/*`, the internal `fcoverage` configuration boundary, downstream scaling-factor loading, README smoothing usage, and existing fragment-count-weight tests in `tests/test_normalize_genome_command.rs`, `tests/test_cli_smoke.rs`, and `src/commands/coverage_weights/coverage_weights_tests.rs`. I did not run tests.

Shared findings that affect this command:

- G-001 in `00_shared_package_notes.md`: shared fragment length ranges can be inverted without a direct error.
- G-004 in `00_shared_package_notes.md`: `--gc-file` lacks shared fail-fast validation for the required `--ref-2bit`.
- G-010 in `00_shared_package_notes.md`: GC correction packages cannot identify the sample or inputs they were built from.
- G-011 in `00_shared_package_notes.md`: scaling-factor TSV metadata is too thin for safe reuse.
- G-012 in `00_shared_package_notes.md`: short final stride bins are length-weighted only in the numerator.
- G-013 in `00_shared_package_notes.md`: smoothing-weight docs claim the inverted scaling factors have mean 1.0.
- G-014 in `00_shared_package_notes.md`: smoothing-weight TSV writes do not explicitly flush the final writer.
- G-015 in `00_shared_package_notes.md`: all-zero smoothing runs fail with a misleading normalization error.
- G-016 in `00_shared_package_notes.md`: smoothing-weight commands cannot match `fcoverage --ignore-gap` segmentation.

## Findings

### FCW-001 - High - The standalone command feature does not enable the shared module it imports

`cmd_fragment_count_weights` depends on `cmd_fcoverage`, but not on `cmd_coverage_weights` ([Cargo.toml](../../Cargo.toml#L36-L39)). The `coverage_weights` module is compiled only when `cmd_coverage_weights` is enabled ([mod.rs](../../src/commands/mod.rs#L7-L8)). `fragment-count-weights` imports both its shared config type and run helper through that gated module ([config.rs](../../src/commands/fragment_count_weights/config.rs#L1-L3), [fragment_count_weights.rs](../../src/commands/fragment_count_weights/fragment_count_weights.rs#L1-L8)).

Impact: default builds include both features, so this can stay hidden. A feature-minimal build such as `--no-default-features --features cli,cmd_fragment_count_weights` should fail because the fragment-count command points at a module its feature does not enable. That makes the release feature matrix unreliable for this command.

Recommended fix:

- If it is acceptable for the fragment-count feature to also expose `coverage-weights`, make `cmd_fragment_count_weights = ["cmd_coverage_weights"]`. Otherwise, move the shared smoothing-weight config/engine into a module gated by `any(feature = "cmd_coverage_weights", feature = "cmd_fragment_count_weights")`.
- If the module is moved, keep the user-facing command modules thin and command-named; the current implementation shape is otherwise reasonable.
- Add a CI feature-matrix check for each released command feature, including `cli,cmd_fragment_count_weights` without default features.

### FCW-002 - Low - Fragment-count outputs still use coverage terminology

The command documentation correctly frames the raw values as fragment-count density/support rather than ordinary coverage ([config.rs](../../src/commands/fragment_count_weights/config.rs#L32-L39), [config.rs](../../src/commands/fragment_count_weights/config.rs#L57-L60)). The shared writer still logs "overlapping position-coverage" and writes columns named `average_pos_coverage` and `average_overlapping_pos_coverage` for the fragment-count output ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L139-L142), [coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L159-L163)). The CLI smoke test also pins those coverage-named columns for `fragment-count-weights` ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L366-L375)).

Impact: the file name says `fragment_counts`, but the primary numeric columns read like coverage-weight output. Since downstream loaders use only `scaling_factor`, this is mostly a human-facing schema clarity issue, but it makes it easier to miss whether a TSV was built from coverage support or fragment-count support.

Recommended fix:

- Rename the shared columns to neutral names like `average_pos_support` and `average_overlapping_pos_support`, or emit command-specific column names while keeping `scaling_factor` stable.
- Add source metadata as part of G-011 so downstream tools do not need to infer the smoothing source from filenames or column wording.
- Update the CLI smoke and parsing tests to pin the intended fragment-count schema.

## Existing Coverage Notes

The command already has a basic integration check that it writes one stride row per chromosome span ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L260-L290)), an identity-smoothing comparison proving fragment-count mode differs from coverage mode for mixed fragment lengths ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L343-L480)), a CLI smoke test for the output path, metadata line, header, and row count ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L319-L383)), and a helper test proving the shared builder uses `LengthNormalizationMode::UnitMass` for fragment-count weights ([coverage_weights_tests.rs](../../src/commands/coverage_weights/coverage_weights_tests.rs#L395-L416)).

The important missing coverage from this review is standalone feature compilation for `cmd_fragment_count_weights`, fragment-count expected values with real smoothing (`bin_size > stride`), final TSV flushing behavior from G-014, all-zero support diagnostics from G-015, `--ignore-gap` parity from G-016, short-final-bin edge behavior from G-012, richer scaling metadata from G-011, and the intended fragment-count output terminology.
