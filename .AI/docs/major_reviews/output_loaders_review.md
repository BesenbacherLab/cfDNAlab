# Rust output loaders review

Date: 2026-06-17

Scope: implemented Rust output loaders under `src/output_loaders/*`, their direct tests in `tests/test_output_loaders_*.rs` and `src/output_loaders/common_tests.rs`, downstream Rust loader tests under `downstream_tests/rust`, and relevant user/internal docs that describe the implemented loader-facing output schemas. I intentionally did not use the older Rust loader plan as authority. I did not execute tests.

## Summary

This file now tracks only unresolved review points or deliberately deferred design notes. Fixed items were removed so the remaining list is easier to scan.

Remaining themes:

- `fcoverage` still needs stricter support-column and field-specific numeric validation if the loader is meant to enforce the current writer contract.
- `fcoverage` provenance is not fully represented by the loaded table header. This is not a loading bug, but it should be documented and optionally enriched if the writer later adds settings metadata.
- Empty-selection behavior and duplicate-name behavior are inconsistent across loaders.
- Zarr label validation and public API governance need tightening before the loaders are treated as stable external Rust APIs.

## Findings

### OL-002 - P1/P2 validation decision - `fcoverage` accepts rows outside the current writer support contract

This finding is about the implemented writer contract, not a biological universal. I am not claiming that eligible bases can never be smaller than `span_positions - blacklisted_positions` in all future designs. They could be smaller if a future fcoverage implementation subtracts reference-N bases, mappability masks, or other base-level masks from the denominator.

What the current code supports is narrower:

- The current fcoverage spec describes summary support as `eligible_positions = span_positions - blacklisted_positions` ([fcoverage_spec.md](../specs/fcoverage_spec.md)).
- The implemented command paths I re-read use blacklist support as the denominator adjustment for these aggregate rows. I did not find implemented evidence that reference-N bases are subtracted from `eligible_positions`.
- The loader currently checks only upper bounds, not the full relationship between support columns.

The loader validates only partial support invariants:

- Windowed summary rows check interval length equals `span_positions`, `blacklisted_positions <= span`, and `eligible_positions <= span`.
- Grouped rows check `blacklisted_positions <= span_positions` and `eligible_positions <= span_positions`.
- Summary stats parse `nonzero_positions` and `covered_fraction` without checking `nonzero_positions <= eligible_positions` or that finite `covered_fraction` is in `[0, 1]`.

Concrete malformed examples currently pass the loader's support checks:

```text
chromosome start end span_positions blacklisted_positions eligible_positions ...
chr1       0     100 100            20                    90                 ...
```

Under the current writer contract, `eligible_positions` should be `80`, not `90`. A grouped row with `span_positions=100`, `blacklisted_positions=20`, and `eligible_positions=90` has the same problem.

Why not just trust the written file: loaders are a trust boundary for external users. Files can be hand-edited, truncated, mixed across versions, or produced by third-party code. If denominator metadata is impossible under the current schema, reporting that at load time is better than silently passing it into downstream analyses.

Recommended fix:

- Enforce `eligible_positions == span_positions - blacklisted_positions` for both windowed and grouped rows where all three fields are present, if this remains the schema contract.
- Enforce `nonzero_positions <= eligible_positions`.
- Enforce finite `covered_fraction` values are in `[0, 1]`. Avoid an exact equality check against `nonzero_positions / eligible_positions` unless writer rounding is accounted for.
- Add malformed-row tests for each invariant.
- If future fcoverage output subtracts N bases or other masks from `eligible_positions`, add an explicit schema/settings field for that mask source and update the loader invariant. Do not silently overload the current columns.

### OL-014 - P3 provenance/documentation gap - `fcoverage` table headers do not encode all command semantics

This is not a data-loading bug. The loaded rows and values can be valid across modes, and the data object should not become a branchy set of separate types just to represent provenance.

The issue is narrower: the table header alone does not encode every command-mode distinction a user may care about.

Examples:

- `Average` and `AverageOnUniqueBases` both write `average_<signal>`.
- `Total` and `TotalOnUniqueBases` both write `total_<signal>`.
- `SummaryStats` and `SummaryStatsOnUniqueBases` share summary-stat headers.
- Length-normalized outputs use `fragment_mass`, but the header alone does not distinguish unit-mass normalization from restored-mean normalization.

The current fcoverage command does not appear to write a settings JSON file, so a settings-loader path would require writer-side metadata first. Today, the stable source of these distinctions is mostly the command-generated filename and user context.

Impact: external Rust users can load two files with materially different grouped-row semantics and see the same public value mode and signal metadata. That can make user-facing summaries under-informative, especially if the object is passed around without its original path or run settings.

Recommended fix:

- Keep the primary loaded object simple.
- Document on `FCoverageValueMode` and `FCoverageSignal` that they describe table headers, not full command provenance.
- Optionally expose lightweight provenance when available, such as parsed canonical filename semantics or future settings JSON metadata. Missing provenance should be `Unknown` or `None`, not guessed.
- If stronger provenance is needed, add a writer-side fcoverage settings JSON first, then add an optional loader path that attaches it.

### OL-005 - P2 design consistency gap - Midpoint `positions(&[])` cannot produce an empty selection

The midpoint selector accepts explicit position indices through `positions(&[usize])`. Other explicit empty selections in the loader family generally work when they can be represented as zero-sized matrices or arrays.

For positions, an empty explicit slice reaches count reading and fails because the reader needs a bounding position span for non-contiguous selections. This is not about empty command output arrays being biologically meaningful. The meaningful case is selection algebra: if user code computes a filter and gets zero positions, returning a `(groups, lengths, 0)` array can be easier to compose than forcing every caller to special-case empty filters.

Impact: callers cannot use a natural filter pipeline where a computed position selection may be empty. The empty-result case is representable as a `(groups, lengths, 0)` array, so the current failure is an implementation/policy choice rather than a storage-schema limitation.

Recommended fix:

- Decide the policy first: either explicit empty selections are valid across loaders, or they are rejected consistently at selector boundaries.
- If empty selections are valid, short-circuit midpoint count selection when any selected axis is empty and return a correctly shaped zero-value `DenseArray3` without touching Zarr.
- If empty selections are invalid, reject them early with a consistent selector error and do not present zero-sized selections as supported elsewhere.
- Add tests for whichever policy is chosen.

### OL-012 - P3 public API consistency - Loader errors expose `anyhow::Result`

The crate defines a public `cfdnalab::Error` and `cfdnalab::Result` for public helpers. The output loaders instead expose `anyhow::Result` throughout their public methods, including `load_lengths_output`, `load_fcoverage_output`, `load_ends_output`, and `load_midpoints_output`.

Impact: this may be acceptable for IO-heavy loaders, but it is a public API commitment. External users cannot pattern-match loader failures the way they can for interval errors. They also inherit `anyhow` as the error surface.

Recommended fix:

- If loaders are intended as a stable public API, either introduce typed loader errors or explicitly document that loader failures are contextual `anyhow` errors.
- Do not do this for backwards compatibility unless that concern is explicitly in scope.

### OL-013 - P3 selection consistency - Explicit empty global-row selections keep `Global` metadata

The common row resolver allows explicit empty row selections because it returns the provided slice without checking non-empty. For global length and end-motif outputs, selecting `rows(&[])` can produce zero count rows while the selected row metadata remains `Global`.

Impact: this is not a data corruption bug, but it makes the metadata less self-describing: a `Global` row metadata variant normally implies one global row, while the selected matrix has zero rows.

Recommended fix:

- Decide whether explicit empty selections are supported across all loaders.
- If yes, document zero-row metadata behavior and add tests.
- If no, reject empty explicit row selections consistently.

### OL-015 - P2/P3 validation gap - Zarr label constraints are not mirrored by Zarr loaders

The Zarr writers reject public labels with control characters because those labels later move into TSVs, data frames, plots, and command-line examples. The end-motif writer applies that rule to concrete motifs, motif groups, group names, and chromosome labels, and additionally requires concrete `motif_ascii` labels to be ASCII. The midpoint writer applies the same label rule to group names.

The loaders do less:

- `ends` concrete motif labels decode `motif_ascii` with `String::from_utf8()` but never check that the bytes are ASCII.
- `ends` and `midpoints` JSON-label readers check `label_field`, label count, and JSON string type, but do not reject control characters.

Impact: malformed Zarr stores can expose labels that the command writers deliberately refuse to create. That becomes an external-user problem once these labels are copied into data frames, TSV exports, plots, or selector strings.

Recommended fix:

- Move a small public-loader-safe label validator into shared code and call it from both Zarr loaders.
- For concrete end motifs, validate `motif_ascii` bytes are ASCII before UTF-8 conversion.
- Add loader tests for non-ASCII concrete motif bytes and control-character JSON labels.

### OL-016 - P3 API consistency - Group-name duplicate policy differs across loaders

`midpoints` rejects duplicate group names during load by building a unique name-to-index map. `lengths` and `ends` load duplicate group names and report the duplicate only when the duplicated name is looked up. The lengths test suite deliberately pins that behavior by loading two `alpha` rows successfully and then expecting `group_index("alpha")` to fail.

This is not automatically wrong, because package loader policy can reasonably allow a duplicated unrelated name without poisoning lookup for a unique name. But the Rust loaders are inconsistent with each other, and the current Rust docs do not state which policy is intentional.

Impact: users get different failure timing and different tolerance for malformed grouped metadata depending on which loader they call. That matters for libraries that want one generic "load grouped output and then select by name" workflow.

Recommended fix:

- Decide whether Rust loaders should be strict command-output validators or permissive name-lookup APIs.
- If strict validation is desired, reject duplicate group names at load for lengths and ends.
- If permissive lookup is desired, consider whether midpoint load-time rejection is too strict, or document it as a midpoint-specific choice.

### OL-017 - P3 docs/API governance - The Rust public API spec omits `output_loaders`

`src/lib.rs` publicly exposes `cfdnalab::output_loaders`, and `src/output_loaders/mod.rs` contains rustdoc that presents the loader entry points as public APIs. The internal Rust public API spec lists stable public categories but does not include `cfdnalab::output_loaders`.

There is also a small rustdoc feature-gating mismatch: the module-level docs list all four loader entry points unconditionally, while the actual modules and re-exports are feature-gated behind `cmd_lengths`, `cmd_fcoverage`, `cmd_ends`, and `cmd_midpoints`.

Impact: contributors do not have a single current policy for whether loader APIs are stable, command-feature-gated public surface, or an experimental convenience module. External users can discover the module in rustdoc without seeing the feature story up front.

Recommended fix:

- Add `cfdnalab::output_loaders` to `.AI/docs/specs/rust_public_api.md`, including the feature-gating policy and whether `anyhow` errors are part of the intentional surface.
- Gate or phrase module rustdoc so no-default-feature users are not shown unavailable entry points as if they are always present.

### OL-018 - P2/P3 validation gap - `fcoverage` numeric value validation is broader than the writer contract

The fcoverage loader uses a generic `parse_f64_field()` that only parses text as `f64`. It does not reject infinities, negative finite values, or field-specific `NaN` values. The command writer has domain-specific behavior: scalar averages with no denominator are `NaN`, scalar totals are raw non-negative sums, and summary-stat zero-support rows keep raw totals at `0.0` while derived statistics become `NaN`.

The loader tests intentionally preserve `NaN` for undefined fcoverage values, so this should not be "reject all NaN". The issue is that the loader has no field-specific rules at all.

Impact: malformed fcoverage TSVs can return `inf`, `-1`, or `NaN` in fields where command output should not have them. Downstream analyses may treat those as real signal rather than reporting a corrupt output file.

Recommended fix:

- Keep the intentional `NaN` cases, but encode them explicitly by value mode and support metadata.
- Reject infinities for all fcoverage numeric value fields.
- Reject negative finite values for coverage/fragment-mass totals, averages, variance, SD, and raw moments.
- Decide whether summary raw totals should be strict writer-compatible (`0.0` when eligible support is zero) or tolerant. Document the choice and test it.

<!--
### OL-010 - Deferred external-docs gap - Public Rust loaders are barely discoverable

This is intentionally deferred for now. The Rust loader surface is public in code, but the external documentation should wait until the API is stable enough to document without creating churn. When that changes, add a short README/website section with feature flags, one minimal example per loader, and current limitations such as fcoverage positional outputs and grouped group-index file handling.
-->

## Coverage and design notes

Existing Rust loader tests cover a lot of parser and selector behavior:

- `lengths`: global/window/group schemas, gzip/zstd, selectors, duplicate/missing names, malformed headers, header-only files, invalid intervals/fractions/column counts, invalid count values, axis contiguity/bounds, and bedGraph-like input rejection.
- `fcoverage`: scalar and summary aggregate schemas, NaN preservation, duplicate selectors, wrong row/value modes, invalid positions, positional-path rejection, header-only files, group-index file loading, and several malformed headers/metadata cases.
- `midpoints`: Zarr metadata, full and selected count reads, duplicate selectors, missing ranges, selector conflicts, schema-version errors, rank errors, non-zero-based axes, axis contiguity/bounds, finite count validation, and duplicate group names.
- `ends`: dense and sparse stores, global/window/group metadata, motif and motif-group axes, selectors, sparse zero-motif output, wrong schema versions, dense shape mismatch, malformed dense/sparse counts, and unsorted sparse coordinates.

The main remaining test weakness is not the number of hand-authored tests. It is that some loaders, especially fcoverage, lack writer-to-loader compatibility tests based on command-generated artifacts. The R/Python fcoverage package-loader work is tracked as a TODO in `.AI/docs/future/fcoverage_package_loader_plan.md`.

## Checked and not currently findings

- I am not flagging fcoverage `NaN` parsing in general. The tests intentionally preserve `NaN` for undefined averages and summary statistics.
- I am not flagging header-only rejection as a bug. The implemented commands generally avoid producing meaningful zero-row outputs, and the loaders consistently reject missing data rows where a command output is expected.
- I am not flagging noncontiguous fcoverage `group_idx` values as malformed by itself. The lower-level grouped writer tests use noncontiguous group indices, so the loader should not assume `group_idx == row_index`.
- I am not flagging `.zarr` suffix requirements. The user accepted that loader paths are cfDNAlab command output paths, not arbitrary renamed Zarr stores.
- I am not using the older loader plan as evidence. All findings above are based on implemented code, implemented tests, or current docs.

## Suggested fix order

1. Tighten fcoverage invariants under the current writer contract (OL-002 and OL-018), including field-specific numeric validation.
2. Decide the empty-selection policy across loaders (OL-005 and OL-013).
3. Tighten Zarr label validation in the Zarr loaders (OL-015).
4. Decide cross-loader group-name duplicate policy (OL-016).
5. Decide and document fcoverage provenance handling only if users need it through the Rust loader (OL-014). If needed, add writer-side settings metadata first.
6. Clarify the Rust public API policy for `output_loaders` (OL-012 and OL-017).
