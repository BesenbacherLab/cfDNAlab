# Command release audit spec

Date: 2026-03-03

## Scope

This audit covers the nine commands documented as released in `README.md`:

- `cfdna fcoverage`
- `cfdna midpoints`
- `cfdna lengths`
- `cfdna gc-bias`
- `cfdna ref-gc-bias`
- `cfdna coverage-weights`
- `cfdna bam-to-bam`
- `cfdna bam-to-frag`
- `cfdna frag-to-bam`

Baseline compile check status:

- `cargo check --features cli,plotters` passes.

## Summary of easy wins

Priority order is based on release risk versus implementation effort.

### P0: user-facing correctness and trust

1. Remove public TODO placeholders in `README.md`.
- Current TODOs are visible in FAQ, recipes, and output format docs.
- This is low effort and high trust impact.

### P1: command-level validation gaps

4. Add command-level integration tests for commands that still have near-none.
- `gc-bias`: internal logic is tested and CLI smoke exists, but deep command-level E2E coverage is still missing.
- `ref-gc-bias`: helper logic is tested and CLI smoke exists, but deep command-level E2E coverage is still missing.

5. Expand smoke-only command tests to cover key behavior toggles.
- `bam-to-frag`: currently mostly one smoke path.
- `midpoints`: currently one integration path.
- `fcoverage`: command-level coverage is now broad. Remaining work is mostly edge-combination saturation rather than missing core behavior tests.

### P2: release process quality

6. Add a command contract table in docs.
- For each command: output files, schema, and minimal expected invariants.
- This prevents drift between CLI behavior and README claims.
Answer: The README is for the external user only. It is a communication document not a AI spec.

1. Add a manual release validation checklist file with pass/fail boxes.
- You already require manual validation.
- Turning this into a strict checklist reduces missed steps.
TODO: Please create an initial plan with such steps

## Command-by-command findings and proposed validation

## 1) cfdna lengths

Current state:

- Strongest test coverage among released commands.
- Includes multiple integration-style checks and edge cases.

Easy wins:

- CLI-level smoke coverage is in place (`--help` coverage and a minimal invocation).

Manual validation targets:

- Verify global output shape and index-to-length mapping.
- Verify windowed output consistency against global counts on a small fixture.
- Verify behavior for `--indel-mode ignore|adjust|skip` on known indel fixtures.

## 2) cfdna bam-to-bam

Current state:

- Good coverage for mapq, blacklist interactions, window overlap behavior, chromosome sort toggle, and scaling tag write.

Easy wins:

- Add unpaired mode coverage (`--reads-are-fragments`).
- Add GC-file validation path tests (missing `--ref-2bit`, invalid package, drop-invalid behavior).
TODO: Since these are similar across commands, we should probably have a standard set of tests that all commands must be tested on?

Manual validation targets:

- Verify FLEN, GC, and COV tags are present and numerically sane when enabled. (should look at distributions and summary stats for this!)
- Verify read count changes across blacklist and window filters.

## 3) cfdna bam-to-frag

Current state:

- One broad smoke integration test exists.

Easy wins:

- Add focused tests for:
  - `--by-bed` filtering behavior.
  - blacklist strategy variants (`any`, `all`, `midpoint`, `proportion`).
  - optional output columns when GC and/or scaling are enabled.
  - unpaired mode.

Manual validation targets:

- Confirm column order and header file match for all correction combinations.
- Confirm output sortedness by `(chromosome, start, end)`.

## 4) cfdna frag-to-bam

Current state:

- Dedicated command-level test suite now covers:
  - parsing failures and invalid strand handling.
  - ordering constraints and chromosome bounds checks.
  - filter behavior (mapq, length, blacklist) and empty-result behavior.
  - BAM header/tid ordering and `bam-to-frag -> frag-to-bam` round-trip restoration.
  - extra-column mapping for `gc_weight`, `scaling_weight`, `flen`.
  - header-source behavior (inline, explicit file, companion file), conflict detection, and unknown-column modes (`--ignore-extras`, `--allow-unknown-extras`).

Easy wins:

- Add explicit stderr warning assertions for `--allow-unknown-extras`.
- Add compressed companion-header detection edge-case tests (`.frag.tsv.gz`, `.frag.tsv.zst`, `.frag.tsv.bgz`).

Manual validation targets:

- Round-trip fixture: `bam-to-frag` -> `frag-to-bam`, then inspect BAM for expected records and header.

## 5) cfdna fcoverage

Current state:

- CLI smoke covers a minimal successful invocation.

- Command tests now cover:
  - basic per-position output.
  - `--keep-zero-runs` toggling for positional output.
  - `--ignore-gap` for paired reads, plus rejection in unpaired mode.
  - unpaired mode basic behavior, plus rejection of `--require-proper-pair`.
  - conservation of total covered bases between positional output and `--by-size total`.
  - tile-size invariance of those totals across multiple tile sizes.
  - `--by-size total` and `--by-size average`, including a non-aligned tile-size average case.
  - `--by-bed average` and `--by-bed total`.
  - `--per-window unique-positions` with overlapping BED windows.
  - `--per-window indexed-positions` with preserved original window indices.
  - multi-run positional output within a single window for both positional BED modes.
  - explicit rejection of `--by-size` with positional per-window modes.
  - blacklist masking in positional output, including a blacklist interval that crosses a tile boundary.
  - blacklist-aware totals in aggregated output, including tile-size invariance when the blacklist crosses a tile boundary.
  - non-integer positional output with scaling and rounding.
  - GC-tag weighting in both unpaired mode and paired-end averaging mode.
  - paired GC-tag edge cases, including invalid-mate fallback and zero-mate zeroing behavior.
  - GC-tag missing or invalid fallback behavior, plus `--drop-invalid-gc`.
  - GC-file validation for missing `--ref-2bit`.
  - GC-file weighting from a tiny reference package, plus `--drop-invalid-gc` fallback/drop behavior.
  - blacklist-aware averages with an explicit masked-denominator check.
  - cross-tile `--by-bed` aggregate reduction with exact output invariance across tile sizes.

- Internal coverage and window result helpers also have dedicated unit tests in `tests/test_coverage.rs`,
  including blacklist-aware positional windows and aggregate calculations.

Easy wins:

- No obvious high-impact easy wins remain in the core `fcoverage` command paths.

- Remaining work is mostly combinatorial saturation if we want even more release confidence, for example:
  - pairing GC-file correction with more aggregate/window modes.
  - checking summary statistics text for GC-failure reporting.
  - adding another blacklist + `--by-bed average` combination test if we want one more denominator check through the BED reducer.

Manual validation targets:

- Spot-check one GC-corrected run and confirm the summary statistics look sane.

- Spot-check one positional output with a very small tile size to make sure the documented segment splitting is understandable.

## 6) cfdna midpoints

Current state:

- One integration test validates output exists, shape, and group mapping.
- Internal accumulator has unit tests.
- Length-bin validation and `start:end:step` parsing support are now implemented and wired into command flow.

Easy wins:

- Add tests for malformed interval files (unsorted, mixed width, empty).
- Add tests for blacklist/scaling/GC interaction paths.

Manual validation targets:

- Confirm group-level profile shape is `(groups, length_bins-1, window_size)`.
- Confirm known interval groups aggregate as expected on synthetic fixtures.

## 7) cfdna coverage-weights

Current state:

- Command tests now cover:
  - output file creation.
  - invalid stride checks.
  - contiguity and chromosome endpoint invariants.
  - deterministic scaling regression values on a simple fixture.

Easy wins:

- Add blacklist and unpaired-mode behavior tests.

Manual validation targets:

- Confirm scaling factor mean approximately 1.0 among non-zero bins.
- Confirm uncovered bins are 0 and covered bins are positive.

## 8) cfdna gc-bias

Current state:

- Strong unit coverage for binning, smoothing, windows, interpolation, and outlier logic.
- CLI-level smoke command coverage exists, but deep command-level E2E tests are still missing.

Easy wins:

- Add end-to-end command test producing package files from tiny BAM + tiny reference assets.
- Validate output package schema and metadata consistency.
- Add failure tests for bad reference package directory contents and invalid outlier args combinations.

Manual validation targets:

- Inspect generated correction matrix for finite values and expected bounds.
- Verify compatibility with downstream commands (`lengths`, `fcoverage`, `midpoints`).

## 9) cfdna ref-gc-bias

Current state:

- Tests cover key counting helpers and CLI-level smoke command path, but deep command-level E2E tests are still missing.

Easy wins:

- Add command-level integration test asserting output artifacts and metadata coherence.
- Add failure tests for invalid smoothing and impossible effective length settings.

Manual validation targets:

- Verify reproducibility with fixed `--seed`.
- Verify expected output changes when changing `--end-offset`.

## Release-level validation matrix

Run these as a fixed gate before release tagging:

1. Build gate
- `cargo check --features cli,plotters`
- `cargo build --release --features cli,plotters --bin cfdna`

2. Help/docs gate
- `cfdna --help`
- `cfdna <command> --help` for all nine release commands.
- No TODO placeholders in README or command help text.

3. Test gate
- Full `cargo test` with default release command feature set.
- Dedicated `frag-to-bam` command tests are green.
- CLI smoke tests for `gc-bias` and `ref-gc-bias` are green.
- Deep command-level E2E tests for `gc-bias` and `ref-gc-bias` remain pending.

4. Manual gate
- Run one small manual invocation per command with fixed tiny fixtures.
- Verify output file names, schema columns, and at least one numeric invariant per command.

## Suggested execution phases

Phase 1 (quick confidence, low risk):

- README TODO cleanup.
- CLI help text corrections.

Phase 2 (highest risk reduction):

- Add deep command-level E2E tests (`gc-bias`, `ref-gc-bias`).
- Expand smoke-only command tests (`bam-to-frag`, `midpoints`, `fcoverage`, `coverage-weights`).

Phase 3 (release hardening):

- Add command output contract table to docs.
- Create and use a strict manual validation checklist for signoff.
