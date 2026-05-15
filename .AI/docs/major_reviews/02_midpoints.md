# Midpoints Command Review

Current review date: 2026-05-14

This is a fresh review of the current `cfdna midpoints` implementation. Older
notes that were fixed by the recent midpoint work have been removed so this file
only tracks still-relevant release risks, design issues, and coverage gaps.

I reviewed:

- `src/commands/midpoints/`
- `src/shared/bed.rs`, `src/shared/interval.rs`, `src/shared/overlaps.rs`,
  `src/shared/window_fetch.rs`, `src/shared/tiled_run.rs`
- `src/shared/blacklist/`, `src/shared/length_axis.rs`,
  `src/shared/scale_genome.rs`
- midpoint-related integration tests in `tests/test_profile_groups_command.rs`,
  `tests/test_tiling.rs`, `tests/test_cli_smoke.rs`, and
  `tests/test_cross_command_artifact_matrix.rs`
- public docs in `README.md` and `.AI/docs/specs/midpoints_spec.md`

I did not run tests, per project instruction.

## Release Triage

No release-blocking findings remain in this review.

Important but not necessarily release-blocking:

- M-001: QC plot x-axes are not the same coordinate contract as the output
  tensor, especially after `--bin-size > 1`.
- M-002: File-based GC correction still builds GC prefixes before no-window tile
  pruning.
- M-003: Midpoint internals are publicly exported as crate API.
- M-004: Several high-value edge cases remain untested.

## Findings

### M-001 - Medium - QC plot x-axes do not match the output coordinate contract

Evidence:

- `src/commands/midpoints/plotting.rs:63-78` builds plot x-values from the final
  number of columns, not from the original interval positions or final
  position-bin widths.
- For odd `window_size`, the plot uses centered coordinates; for even
  `window_size`, it uses raw `0..window_size`.
- `src/commands/midpoints/plotting.rs:100-137` labels both line plots and
  heatmaps as `"Position"`.
- `src/commands/midpoints/postprocess.rs:56-69` allows final binning after
  smoothing and flank trimming, so plot columns can represent position bins
  rather than single bp positions.

Impact:

For `--bin-size > 1`, a plot x-coordinate of `5` means final output column 5,
not genomic offset 5 bp. For odd unbinned windows, the plot is centered even
though the output sidecar says `interval_relative_zero_based`. For even windows,
the plot switches to zero-based. The visual QC can therefore disagree with the
`.npy` coordinate contract.

Recommendation:

- Build a shared position-axis helper from `ProfileLayout` and settings
  semantics.
- Use output-bin intervals or centers in bp units.
- Keep the same convention for odd and even intervals.
- Record or display enough plot metadata that users can distinguish zero-based
  position bins from centered bp offsets.
- Add unit tests for plot-axis construction. The current plotting path has no
  direct tests.

### M-002 - Low - GC reference prefixes are built before no-window tile pruning

Evidence:

- `src/commands/midpoints/midpoints.rs:598-610` reads the 2bit sequence and
  builds GC prefixes for the tile fetch span when GC correction is enabled.
- `src/commands/midpoints/midpoints.rs:632-644` only afterward determines that
  no windows overlap the tile core and returns early.
- `.AI/docs/specs/midpoints_spec.md:80` already records this as an open note.

Impact:

Sparse targeted runs with GC correction can waste reference IO and prefix-build
CPU for tiles that will not count any midpoint sites. This should not change
counts, but it is an avoidable performance cost on sparse interval sets.

Recommendation:

- Move GC-prefix construction after `get_overlapping_sites_and_adapt_fetch_to_extremes`.
- If reading the narrowed `fetch_span`, make sure the GC prefix coordinate origin
  changes with it; `fetch_start` currently assumes `tile.fetch_start()`.

### M-003 - Low - Midpoint internals are exposed as public modules

Evidence:

- `src/commands/midpoints/mod.rs:1-11` exports `counting_by_group`,
  `postprocess`, `smoothing`, and `windows` as `pub mod`.
- `src/commands/midpoints/counting_by_group.rs:201-210` exposes
  `ProfileGroupsCounts` and its fields publicly.

Impact:

This makes internal tensor storage, sparse merge helpers, and postprocessing
types look like supported public API. That raises the compatibility burden for
implementation details that are still changing around release.

Recommendation:

- Make internal modules `pub(crate)` unless they are intentionally part of the
  library API.
- Keep the public surface to configuration, command execution, and any explicit
  reusable settings types.
- If downstream tests need access, prefer crate tests or narrow helper exports
  over making the full internals public.

### M-004 - Low - Remaining high-value coverage gaps

The test suite is much stronger than the older review file suggested. These
items still look worth adding before or soon after release:

- Sparse merge wrap-around: `src/commands/midpoints/counting_by_group.rs:688-695`
  intentionally merges from a non-zero start chunk to the end and then wraps.
  Existing tests cover summing across chunks and valid start-chunk selection, but
  not a forced non-zero wrap-around over entries before and after the split.
- Multiple blacklist inputs: add command tests for more than one blacklist file
  and `--blacklist-min-size`.
- `--require-proper-pair` plus `--reads-are-fragments`: the guard exists in
  `src/commands/midpoints/midpoints.rs:96-98`, but I did not find a midpoint
  command regression test. Other commands test this conflict.
- Non-uniform interval failure through `run`: `ensure_uniform_window_len` exists
  in `src/commands/midpoints/windows.rs:10-38`, but a full command test would
  protect the CLI-facing error path.
- Exact tile-core boundary: add a test where a midpoint falls exactly on
  `core_start` or `core_end` to lock the half-open ownership contract in
  `src/commands/midpoints/midpoints.rs:732-735`.
- Invalid GC neutralization: I did not find a midpoint command test with
  `neutralize_invalid_gc: true`. Other commands cover this behavior.
- Plotting helper tests: no direct tests currently lock plot-axis construction
  or invalid plot-group behavior.

## Stale Notes Removed

The following older review items are no longer active based on the current code
and tests:

- Exact midpoint tensor axis and full expected-count fixture coverage now exists.
- Single-thread versus multi-thread midpoint output equivalence is now tested.
- Group index output includes `eligible_intervals` and has focused tests,
  including zero-count retained groups.
- Smoothing/binning behavior, last-bin width, settings JSON, and default
  `savgol=165` parsing now have direct tests.
- The README midpoint example now mentions default profile smoothing and shows
  `--smoothing none`.
- The README interval description now mentions optional score and strand columns.
- The within-group duplicate/overlap finding was removed. The `--intervals` help
  now states that duplicate intervals and within-group overlaps are user-owned
  clean-ups.
- The settings-sidecar strand-orientation finding was removed. The sidecar is
  intentionally minimal and does not record command invariants.
- The plot-failure finding was removed. The command intentionally moves core
  outputs into place before plotting and still returns an error if QC plotting
  fails.
- Strand-aware command behavior is implemented and tested for forward/reverse
  profile mirroring.
- The current spec records the dense public output and sparse internal tile
  partials.
- The dense allocation/RAM note was removed. The existing size warning is about
  output file size, and dense output naturally requires enough memory for the
  dense merge.
- The blacklist-strategy note was removed. `--blacklist-strategy` is explicitly
  fragment-level, while interval prefiltering is a separate edge-bias filter.
- The strand sampling note was removed. The current heuristic is accepted for
  this command's expected BED-like inputs.
