# `cfdna fcoverage` review

Date: 2026-04-24

Scope: `src/commands/fcoverage/*`, the CLI dispatch for `fcoverage`, and directly used shared helpers. I also read the existing `tests/test_fcoverage_command.rs` function list and targeted sections to understand coverage shape. I did not run tests.

Shared findings that affect this command:

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release docs/API polish:

- F-006: fix blacklist help for positional outputs.

Post-release performance:

- G-006: sparse-window GC reference pruning.

## Findings

### F-006 - Low - Blacklist help says positional values become `f32::NaN`, but positional writers omit masked bases

The config long help says blacklisted positions are set to `f32::NaN` ([config.rs](../../src/commands/fcoverage/config.rs#L66-L69)). The actual positional writer skips masked bases and splits runs around them ([writers.rs](../../src/commands/fcoverage/writers.rs#L1179-L1241)). For aggregate outputs, the implementation excludes masked bases from sums/averages instead of writing NaNs.

Impact: the user-facing text describes an internal/statistical idea, not the emitted bedGraph/TSV behavior. This matters because bedGraph consumers will see gaps, not rows containing `NaN`.

Recommended fix:

- Change the help text to say: positional outputs omit blacklisted bases; aggregate outputs exclude blacklisted bases from eligible sums/averages and report `blacklisted_positions`.

## Existing coverage notes

The command already has unusually broad end-to-end coverage in `tests/test_fcoverage_command.rs`, including normalization, restore-mean, gap handling, blacklist masking, tile-boundary invariance, by-size/by-BED/grouped-BED modes, GC file/tag paths, and invalid mode combinations. The deferred sparse-window GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.
