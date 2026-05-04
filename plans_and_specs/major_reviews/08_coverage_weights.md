# `cfdna coverage-weights` review

Date: 2026-04-24

Scope: `src/commands/coverage_weights/*`, the internal `fcoverage` configuration boundary, downstream scaling-factor loading in `src/shared/scale_genome.rs`, README smoothing usage, and existing coverage-weight tests in `tests/test_normalize_genome_command.rs`, `tests/test_cli_smoke.rs`, `tests/test_cross_command_artifact_matrix.rs`, and `src/commands/coverage_weights/coverage_weights_tests.rs`. I did not run tests.

Shared findings that affect this command:

- Original review note: no active shared correctness findings were tracked in `00_shared_package_notes.md` at the time. Current re-review additions below list shared findings that now affect this command.

## Release triage

Pre-release correctness/safety:

- None active in the original review pass; see re-review additions below.

## Findings

No additional `coverage-weights`-specific findings were active in the original review pass.

The command's feature wiring is internally consistent in the current feature matrix: `cmd_coverage_weights` enables `cmd_fcoverage` ([Cargo.toml](../../Cargo.toml#L36-L39)), the module is gated by `cmd_coverage_weights` ([mod.rs](../../src/commands/mod.rs#L7-L8)), and the CLI subcommand is gated by the same feature ([cli_app.rs](../../src/cli_app.rs#L5-L6), [cli_app.rs](../../src/cli_app.rs#L79-L80)). The command is also a thin wrapper over the shared smoothing implementation, with coverage mode selected by `normalize_by_length = false` ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L89-L90), [coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L233-L238)).

## Existing Coverage Notes

The command already has direct coverage for output row ranges and non-zero support behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L208-L258)), contiguous output bins ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L293-L341)), hand-derived smoothing and scaling values for a simple fragment ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L484-L563)), unpaired read-as-fragment behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L565-L640)), bin-size validation ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L642-L689)), support-floor and length-weighted global mean behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L691-L907)), one shared global mean across chromosomes ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L909-L1059)), default MAPQ behavior ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L1061-L1274)), chromosome row ordering ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L1276-L1373)), blacklist masking ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L1375-L1457)), CLI smoke output shape ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L253-L317)), and cross-command consumption of a real scaling TSV with GC correction ([test_cross_command_artifact_matrix.rs](../../tests/test_cross_command_artifact_matrix.rs#L96-L147), [test_cross_command_artifact_matrix.rs](../../tests/test_cross_command_artifact_matrix.rs#L240-L450)).

The important missing coverage from this review is standalone feature-matrix compilation.

## Released-command re-review additions (2026-05-04)

### Shared findings that affect this command

- G-019: `coverage-weights` runs internal `fcoverage`, so it inherits the tiled raw-chromosome temporary filename issue from that path.
- G-021: `coverage-weights` exposes shared GC correction options, including `--gc-tag`.
- G-022: the final scaling TSV and the internal `fcoverage` source directory both use the unchecked output prefix.

G-023 can affect the scientific validity of a provided GC correction package indirectly, but the correction-package identity fix belongs upstream in `gc-bias`.

### Release triage additions

Pre-release correctness/safety:

- CW-001: partially blacklisted stride bins are weighted as full-length support during smoothing and global normalization.
- G-021: overlong `--gc-tag` values should fail before reading BAM records.
- G-022: unchecked output prefixes can write files outside the requested output directory.
- G-019: raw chromosome names can escape temporary filename boundaries through the internal `fcoverage` run.

### CW-001 - Medium - Partially blacklisted stride bins are treated as full-length support

`coverage-weights` builds its raw stride values through internal `fcoverage --by-size <stride> --per-window average` ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L42-L54), [coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L272-L292)). In masked mode, `fcoverage` computes the average over eligible positions only, returning `NaN` only when a window has zero eligible positions ([tiling.rs](../../src/commands/fcoverage/tiling.rs#L590-L614)).

That is correct for the raw stride value, but the next stage loses the denominator. The scaling-weight loader parses `blacklisted_positions` and discards it, storing only the average coverage value ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L340-L366)). The triangular smoother then weights each finite stride average by the full stride interval length, not by the number of eligible bases that supported that average ([striding.rs](../../src/commands/coverage_weights/striding.rs#L160-L175)). Global normalization uses the same full interval length for finite non-zero smoothed rows ([striding.rs](../../src/commands/coverage_weights/striding.rs#L212-L239)).

Impact: a stride that is mostly blacklisted but has a small unmasked sliver influences smoothing and the global mean as much as a fully eligible stride of the same genomic span. This can distort scaling factors near partially masked regions. Fully blacklisted strides are already handled because their average is `NaN`; the gap is partial support.

Recommended fix:

- Carry an eligible-base support length in `StrideBin`, derived as `end - start - blacklisted_positions` for coverage mode.
- In `fill_triangular_overlap()` and `normalize_average_overlap_by_global_mean()`, weight finite coverage averages by eligible support length rather than full interval length.
- Add a regression with a mostly blacklisted stride and a fully eligible stride proving the smoothed value and global mean are eligible-base weighted. Keep the existing fully blacklisted regression.
