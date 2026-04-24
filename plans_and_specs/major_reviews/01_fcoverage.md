# `cfdna fcoverage` review

Date: 2026-04-24

Scope: `src/commands/fcoverage/*`, the CLI dispatch for `fcoverage`, and directly used shared helpers. I also read the existing `tests/test_fcoverage_command.rs` function list and targeted sections to understand coverage shape. I did not run tests.

Shared findings that affect this command:

- G-001 in `00_shared_package_notes.md`: inverted fragment length ranges are not rejected.
- G-002 in `00_shared_package_notes.md`: README examples are not clean runnable snippets.
- G-003 in `00_shared_package_notes.md`: tiled temp directories need cleanup guards.
- G-004 in `00_shared_package_notes.md`: `--gc-file` needs shared fail-fast `--ref-2bit` validation.
- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Findings

### F-001 - High - Positional output writers discard write errors

`write_bedgraph_runs()` and `write_windowed_runs()` call `writeln!` inside a callback and assign the result to `_` ([writers.rs](../../src/commands/fcoverage/writers.rs#L981-L1003), [writers.rs](../../src/commands/fcoverage/writers.rs#L1040-L1073)). The comment says errors are "bubbled up by caller on flush", but that is not a reliable contract for every `Write` implementation, and the functions return `Ok(())` even when an individual row write failed.

Affected outputs:

- whole-genome positional bedGraph
- `--by-bed --per-window unique-positions`
- `--by-bed --per-window indexed-positions`
- restore-mean positional paths before final scaled merge

Impact: disk-full, interrupted pipe, or writer failures can be reported as successful output or can produce truncated positional files without the failing write being attached to the row that failed.

Recommended fix:

- Make `visit_runs_in_window()` return `Result<()>`, with a callback returning `Result<()>`.
- Propagate each `writeln!` with `?`.
- Add a tiny writer test with a custom failing writer so this cannot regress without needing filesystem tricks.

### F-002 - High - Covered by shared GC validation and temp cleanup findings

This command is affected by G-003 and G-004 in `00_shared_package_notes.md`. The fcoverage-specific evidence is that `run_inner()` reaches temp directory creation before the tile-level `ref_2bit` check, and its cleanup block runs only on the success tail ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L272-L296), [fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L745-L752), [fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L780-L784)).

### F-003 - Medium - Covered by shared sparse-window GC reference finding

The fcoverage-specific evidence for this is now tracked in G-006 in `00_shared_package_notes.md` to avoid repeating the same reference-preload issue across every tiled command.

### F-004 - Medium - `--normalize-by-length=off` is accepted even though it is not a documented user mode

`LengthNormalizationMode` derives `clap::ValueEnum`, and the internal default variant `Off` is not skipped ([config.rs](../../src/commands/fcoverage/config.rs#L17-L23)). The CLI argument uses that value enum with `default_value_t = LengthNormalizationMode::Off` and `default_missing_value = "unit-mass"` ([config.rs](../../src/commands/fcoverage/config.rs#L140-L152)). The help text documents only `unit-mass` and `restore-mean` as user-facing modes ([config.rs](../../src/commands/fcoverage/config.rs#L124-L136)).

Impact: users can spell an explicit no-op mode that the docs do not describe. That is not a correctness failure, but it adds avoidable CLI surface area while the package is still pre-public.

Recommended fix:

- Mark `Off` as skipped for clap value parsing, while keeping it as the Rust default.
- Add a CLI parser regression test that rejects `--normalize-by-length=off` and accepts bare `--normalize-by-length`.

### F-005 - Medium - Fully masked averages are `0` in average mode but `NaN` in summary-stats mode

The scalar aggregate path intentionally returns `0.0` when a masked average has zero eligible positions ([tiling.rs](../../src/commands/fcoverage/tiling.rs#L575-L613)). The summary-stats path returns `NaN` for average, variance, SD, CV, and covered fraction when `eligible_positions == 0` ([writers.rs](../../src/commands/fcoverage/writers.rs#L153-L178)).

Impact: the same fully blacklisted window can mean "zero coverage" in `--per-window average` but "undefined average" in `--per-window summary-stats`. The `blacklisted_positions` column lets careful users detect this, but downstream code that consumes only the value column can silently treat missing data as true zero.

Recommended fix:

- Decide the public semantic intentionally. I would lean toward `NaN` for scalar average when the denominator is zero, but that may require checking downstream consumers.
- If keeping `0.0`, document the special case directly in the `--per-window average` help text.
- Add a command-level test that pins the chosen behavior for fully blacklisted windows.

### F-006 - Low - Blacklist help says positional values become `f32::NaN`, but positional writers omit masked bases

The config long help says blacklisted positions are set to `f32::NaN` ([config.rs](../../src/commands/fcoverage/config.rs#L66-L69)). The actual positional writer skips masked bases and splits runs around them ([writers.rs](../../src/commands/fcoverage/writers.rs#L1179-L1241)). For aggregate outputs, the implementation excludes masked bases from sums/averages instead of writing NaNs.

Impact: the user-facing text describes an internal/statistical idea, not the emitted bedGraph/TSV behavior. This matters because bedGraph consumers will see gaps, not rows containing `NaN`.

Recommended fix:

- Change the help text to say: positional outputs omit blacklisted bases; aggregate outputs exclude blacklisted bases from eligible sums/averages and report `blacklisted_positions`.

## Existing coverage notes

The command already has unusually broad end-to-end coverage in `tests/test_fcoverage_command.rs`, including normalization, restore-mean, gap handling, blacklist masking, tile-boundary invariance, by-size/by-BED/grouped-BED modes, GC file/tag paths, and invalid mode combinations. The most important fcoverage-specific missing tests from this review are writer-error propagation, CLI parser behavior for `--normalize-by-length=off`, and the chosen scalar-average behavior for fully masked windows. Shared missing coverage for early `--gc-file`/`--ref-2bit` validation and sparse-window GC reference pruning is tracked in `00_shared_package_notes.md`.
