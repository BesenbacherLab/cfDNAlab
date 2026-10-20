# Outliers feature audit

2026-10-20. Scope: Cargo dependencies, command registration, the library API, CLI parsing and
rendering, plotting, output registration, tests, documentation and CI. Comparisons use GC bias,
lengths and midpoints. Current invariants are in [the command spec](../specs/outliers.md#cargo-features).

## Findings and changes

| Area | Finding | Resolution |
| --- | --- | --- |
| Cargo | `cmd_outliers` forced `plotters` on, unlike the other commands. | Remove that dependency. Plotters remains optional and enabled by the existing defaults. |
| Compilation | PNG modules, imports and preserved-context BAM rereads were unconditional. | Compile them only with `plotters`. Keep all numerical analysis and text writers available. |
| Output contract | The manifest and Rust result promised a global PNG in every build. | List the PNG only when produced. `output_plot` is now `Option<PathBuf>`, with `None` without plotting. |
| Explicit requests | A build without plotting had no useful handling for `--plot-preserved-regions`. | Reject the request during configuration validation, before opening the BAM. |
| Reporting | Text outputs did not identify plotting availability. | Add `plotting_available` metadata. Preserve all inspection-table coordinates and use `NA` for unavailable plot files. |
| Test coverage | Command output tests assumed PNGs, and the public rendering/help tests omitted outliers. | Cover manifests, result paths, requests, stale PNGs, CLI presence and rendering with the relevant feature combinations. Keep the numerical expectations unchanged. |
| CI | Per-command checks and the non-CLI test feature list omitted outliers. | Add it to those existing lists, using the same checks as the other commands. |

Old PNGs are left intact when a prefix is reused, consistent with the existing optional-plot
behavior. They are excluded from a new text-only run's manifest, and its metadata identifies
plotting as unavailable. A dedicated regression test covers this case.

## Feature responsibilities

- `cmd_outliers` enables the implementation and the `run_like_cli::outliers` API. Its
  `cmd_fcoverage` dependency is needed for shared fragment coverage construction and counters.
  This transitively enables GC support, matching fcoverage's existing dependency structure.
- `argmin`, `argmin-math` and `statrs` are optional dependencies activated by `cmd_outliers`.
  They support ZINB optimization and probability calculations, independently of PNG output.
- `cli` enables Clap derives, parsers, command registration and dispatch. Config construction
  and equivalent-command rendering work without it. Module and public API registration follow
  `cmd_outliers`, as for the other commands.
- `plotters` enables the global PNG and optional preserved-region PNGs. Plotting controls remain
  visible in configs and help, following midpoints' approach. Outliers explicitly rejects an
  unsupported PNG request instead of silently ignoring it. Context options remain useful without
  Plotters because their coordinates are included in the inspection table.
- `testing` exposes fixture builders for integration tests. `docs_gen` remains independent of
  plotting. Existing internal build-script groups already cover the BAM, tiling, temporary-file
  and progress helpers needed by outliers.

No model-selection defaults, command defaults, counting rules or biological decisions changed.
User-facing guides, CLI help and the generated README assume default features. Non-default
build details remain in this internal audit and the command spec.

## Verification

Compilation checks use `cargo check --offline --locked`. Tests were compiled, never executed.

| Build | Result |
| --- | --- |
| `--no-default-features --features cmd_outliers`, with and without `cli` | Passed without warnings. |
| Outliers and public rendering integration tests, all four CLI/Plotters combinations | Passed without warnings. |
| Unit and integration tests without Plotters, with the broader helper command set, CLI on and off | Passed. Warnings are in fragment-kmer code and existing test helpers. |
| Default build and `--tests --features testing,cmd_outliers` | Passed without warnings. |
| `--all-features` | Passed. Warnings are in fragment-kmer and visualize-positions code. |
| fcoverage-only library and public CLI-rendering control | Passed without warnings. The outliers command is excluded from this feature set. |
| `cmd_fcoverage,plotters` with outliers disabled | Passed without warnings. |
| Library without any command features | Passed with shared-code unused-item warnings. |
| Docs example with `--no-default-features --features cli,docs_gen,cmd_outliers` | Passed without warnings. |
| Rust formatting, diff whitespace and workflow YAML syntax | Passed. |

The minimal integration checks select `test_outliers_command`, `test_cli_command_rendering`
and `test_cli_command_rendering_no_cli`. They include the existing hand-derived numerical
assertions as well as the new feature behavior. The broader unit-test feature set is:

```text
testing,cmd_outliers,cmd_coverage_weights,cmd_fragment_count_weights,cmd_ends,cmd_lengths,
cmd_midpoints,cmd_ref_kmers,cmd_bam_to_bam,cmd_fragment_kmers
```

Raw commands and results are recorded locally in
`/private/tmp/cfdnalab-outliers-features-results.json`, with one log per check.

The dependency trees already confirm that `cmd_outliers` without CLI or plotting has `argmin`,
`argmin-math` and `statrs`, and no Plotters, font backend or Clap dependencies. A fcoverage-only
control has none of those optional dependencies.

## Separate shared-test limitation

`cargo check --no-default-features --tests --features testing,cmd_outliers` still fails because
shared unit tests reference helpers belonging to disabled commands. The fcoverage-only control,
with outliers disabled, produces the same 39 errors in the same shared files. Examples include
`cli_common_tests.rs`, `output_loaders/common_tests.rs`, and the shared overlap, window-fetch,
scaling and BAM tests. Coverage-weights tests also reference the feature-gated fragment-count
variant unconditionally when that companion command is absent.

This is a repository test-gating issue, not a dependency needed by outliers. The product's
minimal builds and its minimal public integration-test builds pass. Outliers' private tests
also compile without Plotters using the broader command set above.

Follow-up: match each affected shared test's cfg to the production helper it exercises, while
retaining valid shared coverage for smaller builds. Do not enable extra product dependencies
solely to make the shared test harness compile. This audit leaves those unrelated test modules
unchanged. CI uses the existing shared command checks.
