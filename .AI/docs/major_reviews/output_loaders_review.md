# Rust output loaders review

Date: 2026-06-17

Last updated: 2026-06-18

Scope: implemented Rust output loaders under `src/output_loaders/*`, their
direct tests in `tests/test_output_loaders_*.rs` and
`src/output_loaders/common_tests.rs`, downstream Rust loader tests under
`downstream_tests/rust`, and relevant user/internal docs that describe the
implemented loader-facing output schemas.

There are no remaining active findings from this review.

## Recently Fixed

- `fcoverage` support metadata now enforces
  `eligible_positions = span_positions - blacklisted_positions` where those
  fields are present.
- `fcoverage` summary rows now reject `nonzero_positions > eligible_positions`,
  finite `covered_fraction` outside `[0, 1]`, infinities, and negative finite
  signal values where the writer contract does not allow them. Documented
  zero-support `NaN` values are preserved.
- `fcoverage` exposes lightweight filename metadata for canonical action and
  length-normalization filename parts, with `Unknown` when a renamed file does
  not carry those parts.
- Midpoint outputs with empty count axes are rejected at load, and explicit
  `positions(&[])` selections now fail at the selector boundary.
- Public loader methods now return `OutputLoaderResult<T>` and
  `OutputLoaderError` instead of exposing `anyhow::Result`.
- Explicit `rows()` selection on global lengths and end-motif outputs now fails
  instead of returning global metadata with a selected row axis.
- Duplicate grouped labels are rejected during load for lengths, ends, and
  midpoints. Grouped fcoverage rows reject duplicate `group_idx`.
- End-motif and midpoint Zarr loaders now reject control characters in public
  JSON labels. End-motif concrete `motif_ascii` labels now require ASCII bytes.
- `cfdnalab::output_loaders` is now listed in the Rust public API spec, and the
  module docs describe the command-feature gates.

## Checked And Not Currently Findings

- I am not flagging fcoverage `NaN` parsing in general. The tests intentionally
  preserve `NaN` for undefined averages and zero-support summary statistics.
- I am not flagging header-only rejection as a bug. The implemented commands
  generally avoid producing meaningful zero-row outputs, and the loaders
  consistently reject missing data rows where a command output is expected.
- I am not flagging noncontiguous fcoverage `group_idx` values as malformed by
  itself. The lower-level grouped writer tests use noncontiguous group indices,
  so the loader should not assume `group_idx == row_index`.
- I am not flagging `.zarr` suffix requirements. Loader paths are cfDNAlab
  command output paths, not arbitrary renamed Zarr stores.
