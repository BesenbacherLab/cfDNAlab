# `cfdna fragment-count-weights` review

Date: 2026-04-24

Scope: `src/commands/fragment_count_weights/*`, the shared smoothing-weight implementation in `src/commands/coverage_weights/*`, the internal `fcoverage` configuration boundary, downstream scaling-factor loading, README smoothing usage, and existing fragment-count-weight tests in `tests/test_normalize_genome_command.rs`, `tests/test_cli_smoke.rs`, and `src/commands/coverage_weights/coverage_weights_tests.rs`. I did not run tests.

Shared findings that affect this command:

- Original review note: no active shared correctness findings were tracked in `00_shared_package_notes.md` at the time. Current re-review additions below list shared findings that now affect this command.

## Release triage

Pre-release correctness/safety:

- None active in the original review pass; see re-review additions below.

## Findings

No active findings in the original review pass.

## Existing Coverage Notes

The command already has a basic integration check that it writes one stride row per chromosome span ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L260-L290)), an identity-smoothing comparison proving fragment-count mode differs from coverage mode for mixed fragment lengths ([test_normalize_genome_command.rs](../../tests/test_normalize_genome_command.rs#L343-L480)), a CLI smoke test for the output path, metadata line, header, and row count ([test_cli_smoke.rs](../../tests/test_cli_smoke.rs#L319-L383)), and a helper test proving the shared builder uses `LengthNormalizationMode::UnitMass` for fragment-count weights ([coverage_weights_tests.rs](../../src/commands/coverage_weights/coverage_weights_tests.rs#L395-L416)).

The important missing coverage from this review is fragment-count expected values with real smoothing (`bin_size > stride`).

## Released-command re-review additions (2026-05-04)

### Shared findings that affect this command

- G-019: `fragment-count-weights` runs internal `fcoverage`, so it inherits the tiled raw-chromosome temporary filename issue from that path.
- G-021: `fragment-count-weights` exposes shared GC correction options, including `--gc-tag`.
- G-022: the final scaling TSV and the internal `fcoverage` source directory both use the unchecked output prefix.

G-023 can affect the scientific validity of a provided GC correction package indirectly, but the correction-package identity fix belongs upstream in `gc-bias`.

### Release triage additions

Pre-release correctness/safety:

- FCW-001: fully blacklisted stride bins are treated as finite zero-mass bins in fragment-count mode, so they can dampen smoothing and receive nonzero scaling factors.
- G-021: overlong `--gc-tag` values should fail before reading BAM records.
- G-022: unchecked output prefixes can write files outside the requested output directory.
- G-019: raw chromosome names can escape temporary filename boundaries through the internal `fcoverage` run.

### FCW-001 - Medium - Fully blacklisted stride bins are treated as real zero-mass support

`fragment-count-weights` delegates to `run_with_fcoverage()` with length normalization enabled and `ScalingWeightsCommand::FragmentCount` ([fragment_count_weights.rs](../../src/commands/fragment_count_weights/fragment_count_weights.rs#L7-L14)). That command variant makes the internal `fcoverage` run write `total_coverage` rows, not `average_coverage` rows ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L42-L54), [coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L272-L292)).

For masked windows, `fcoverage` excludes blacklisted bases when computing the coverage sum and records how many bases were blacklisted ([tiling.rs](../../src/commands/fcoverage/tiling.rs#L381-L421)). But `CoverageWindowAction::Total` finalizes the value as the raw sum even when there are zero allowed positions ([tiling.rs](../../src/commands/fcoverage/tiling.rs#L590-L614)). A fully blacklisted stride therefore enters the internal TSV as `total_coverage = 0`, not `NaN`.

The scaling-weight loader parses `blacklisted_positions` but discards it, storing the parsed `total_coverage` as the stride value ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L340-L366)). The triangular smoother skips only non-finite stride values, so a fully blacklisted stride with value `0` contributes real zero mass and real kernel weight to neighboring smoothed bins ([striding.rs](../../src/commands/coverage_weights/striding.rs#L160-L175)). Later normalization also treats a raw zero stride as usable when smoothing gives it nonzero neighbor support, which can produce a nonzero scaling factor for the blacklisted row itself ([striding.rs](../../src/commands/coverage_weights/striding.rs#L212-L275)).

Impact: when `fragment-count-weights` is run with blacklists, fully blacklisted stride bins can depress nearby smoothed fragment-mass estimates and can receive usable scaling factors despite having no eligible bases. `coverage-weights` is less exposed to this exact path because its internal `average_coverage` output becomes `NaN` for zero-eligible masked windows.

Recommended fix:

- When loading internal `fcoverage` rows for `ScalingWeightsCommand::FragmentCount`, convert rows where `blacklisted_positions == end - start` to `NaN` before smoothing.
- Keep zero-valued but eligible rows as real zero support.
- Add a focused fragment-count regression with one fully blacklisted stride between supported strides, proving the blacklisted row is skipped during smoothing and receives scaling factor `0`.
