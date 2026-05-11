# G-006 Plan: Origin-Aware GC Prefixes and Safe Reference Pruning

Status: deferred post-release performance optimization. This plan is intentionally not part of the first-release correctness queue unless later evidence shows a correctness bug independent of the pruning/refactor work.

## Summary

Fix G-006 by making GC prefix arrays carry an explicit mapping between reference coordinates and prefix-array coordinates. Then move reference reads after window/fetch pruning and build prefixes only for the narrowed span.

The key safety rule is: sequence bytes do not carry coordinates. Any prefix object built from those bytes must be given the reference-coordinate interval and the prefix-array coordinate interval explicitly, then downstream code must ask the prefix object to convert reference intervals into prefix-local intervals.

Scope includes all current tiled commands with this pattern: `fcoverage`, `wps`, `midpoints`, `ends`, `lengths`, `fragment-kmers`, and `ref-gc-bias`.

## Prefix Coordinate API

- Add a small coordinate-map type, for example:

  ```rust
  pub struct GCPrefixCoordinateMap {
      pub reference_interval: Interval<u64>,
      pub prefix_interval: Interval<usize>,
  }
  ```

  Contract:
  - `reference_interval` is the absolute half-open reference interval that the byte slice represents.
  - `prefix_interval` is the half-open coordinate interval used to index the prefix arrays.
  - For this implementation, `prefix_interval` must be exactly `[0, seq.len())`. Do not introduce non-zero prefix-array origins unless all local count APIs are redesigned around that.
  - Both intervals must be supplied by the caller because the byte slice itself does not carry either coordinate system.

- Extend `GCPrefixes`:

  ```rust
  pub struct GCPrefixes {
      pub coordinates: GCPrefixCoordinateMap,
      pub gc: Vec<u32>,
      pub acgt: Vec<u32>,
  }
  ```

- Keep `build_gc_prefixes(seq)` for origin-zero/local use. It should construct:

  ```rust
  reference_interval = [0, seq.len())
  prefix_interval = [0, seq.len())
  ```

- Add a constructor that requires the coordinate map explicitly:

  ```rust
  pub fn build_gc_prefixes_with_coordinates(
      seq: &[u8],
      coordinates: GCPrefixCoordinateMap,
  ) -> Result<GCPrefixes>
  ```

  Validation:
  - `coordinates.reference_interval.len() == seq.len()`
  - `coordinates.prefix_interval == Interval::new(0, seq.len())?`
  - prefix vectors have length `seq.len() + 1`

- Add conversion/count methods with names that describe coordinate input, not “absolute counts”:
  - `prefix_interval_for_reference_interval(reference_interval) -> Result<Option<Interval<usize>>>`
  - `gc_count_for_reference_interval(reference_interval) -> Result<Option<u32>>`
  - `acgt_count_for_reference_interval(reference_interval) -> Result<Option<u32>>`
  - Optional: `gc_acgt_count_for_reference_interval(reference_interval) -> Result<Option<(u32, u32)>>`

- Keep local methods:
  - `gc_count(prefix_interval)`
  - `acgt_count(prefix_interval)`
  - `get_gc_integer_percentage_for_window(prefixes, prefix_interval, ...)`

  Do not remove these because `ref-gc-bias` and low-level tests intentionally use prefix-local coordinates.
  Local methods keep the current contract: invalid prefix-array intervals are errors. Reference-coordinate methods use `Option` because an interval can legitimately lie outside the loaded reference slice.

## Command Changes

- For `fcoverage`, `wps`, `midpoints`, `ends`, `lengths`, and `fragment-kmers`:
  - Compute the narrowed `fetch_span` before any GC-prefix reference read.
  - If no fetch is needed, return the existing empty/no-output tile result before reading `ref_2bit`.
  - Read GC-correction sequence from `fetch_span.start()..fetch_span.end()`.
  - Build prefixes with:

    ```rust
    GCPrefixCoordinateMap {
        reference_interval: fetch_span,
        prefix_interval: Interval::new(0, seq_bytes.len())?,
    }
    ```

  - Replace manual `shift_left(tile.fetch_start())` GC-correction code with reference-coordinate correction methods.
  - Update any fragment in-bounds guard used for GC correction to compare against `fetch_span.start()` and `fetch_span.end()`, not `tile.fetch_start()` and `tile.fetch_end()`.

- Add corrector wrappers:
  - `GCCorrector::correct_fragment_for_reference_interval(fragment_interval, prefixes)`
  - `LengthAgnosticGCCorrector::correct_fragment_for_reference_interval(fragment_interval, prefixes)`

  These should:
  - take fragment intervals in reference coordinates,
  - contract by `end_offset`,
  - ask `GCPrefixes` to convert to prefix-local coordinates,
  - preserve the existing error when `contract(end_offset)` makes the fragment empty,
  - return `Ok(None)` if the contracted GC window is outside the prefix map,
  - otherwise behave exactly like existing local `correct_fragment`.

- For `gc-bias` observed-fragment counting:
  - Replace separate `sequence_interval` arguments with `prefixes.coordinates.reference_interval`.
  - Use `gc_count_for_reference_interval` and `acgt_count_for_reference_interval` where current code manually shifts by `sequence_interval.start()`.
  - Add an equivalence regression for both `get_fragment_gc` and `set_window_acgt_in_observed_interval`, because these are the two places where `sequence_interval` is currently threaded separately.

- For `ref-gc-bias`:
  - Build/prune tile windows before reference sequence loading.
  - If no BED windows survive for the tile, return empty counts before `read_seq_in_range`.
  - Continue passing prefix-local windows/starts to `count_reference_gc_and_length_by_window`.
  - Build prefixes with an explicit `GCPrefixCoordinateMap` so the loaded reference slice is still documented and testable.
  - Keep the existing `seq_start = tile.fetch_start().min(tile.core_start())` and `seq_end = tile.fetch_end().min(chrom_len as u32)` bounds unless a separate proof shows they are wrong. This plan is not a tile-boundary cleanup.

- Whole-chromosome non-tiled users (`bam_to_frag`, `bam_to_bam`) are not part of G-006 pruning, but must compile. They can keep `build_gc_prefixes(read_seq(...))`, which means origin-zero chromosome coordinates.

## Tests

Testing goal: the suite must distinguish the intended coordinate behavior from the two known
wrong implementations:

- reading GC prefixes from the full tile before window/fetch pruning
- reading the narrowed reference slice but indexing it as though prefix coordinates still start at
  `tile.fetch_start()`

Tests that depend on new API names are refactor-gated. Do not mark G-006 complete until those tests
are added together with the API they exercise.

- Prefix coordinate-map tests:
  - Constructor rejects mismatched `reference_interval` length.
  - Constructor rejects any `prefix_interval` other than `[0, seq.len())`.
  - Non-zero `reference_interval` plus `[0, seq.len())` prefix interval converts reference intervals correctly.
  - Reference intervals outside the map return `None`.
  - Local `gc_count` / `acgt_count` still work unchanged.
  - Status: refactor-gated, because `GCPrefixCoordinateMap` and reference-coordinate count methods
    do not exist before this change.

- Corrector tests:
  - `correct_fragment_for_reference_interval` matches existing `correct_fragment` when using equivalent local intervals.
  - Same test for `LengthAgnosticGCCorrector`.
  - Empty contracted GC windows still return the existing error, not `Ok(None)`.
  - Contracted GC windows outside the prefix map return `Ok(None)`.
  - Include a non-zero `end_offset` fixture where the untrimmed fragment has one GC bin and the
    contracted GC window has another, so the test proves the wrapper trims before converting.
  - Status: refactor-gated, because the reference-coordinate corrector wrappers do not exist before
    this change.

- Current coordinate-helper regression tests:
  - `gc-bias::get_fragment_gc` with `sequence_interval = [900, 961)` and a prefix-local sequence
    proves fragment GC uses `sequence_interval` as the origin.
  - `gc-bias::get_fragment_gc` returns `None` when the contracted fragment lies outside the loaded
    sequence interval.
  - `gc-bias::set_window_acgt_in_observed_interval` proves ACGT support uses `sequence_interval` as
    the origin.
  - `gc-bias::set_window_acgt_in_observed_interval` clips observed support to the loaded sequence
    and errors when there is no overlap at all.
  - `ref-gc-bias::process_tile` masks blacklist intervals against a non-zero loaded-reference
    origin and proves the masked reference coordinates, not prefix-local coordinates, are used.

- No-reference-read pruning tests:
  - For each affected command family, create a no-window or empty-candidate tile with file-based GC correction enabled and a reference-read sentinel.
  - Expected result: the tile returns its normal empty result instead of failing while opening/reading reference.
  - Cover separately:
    - `fcoverage`/`wps` shared fetch-helper family,
    - `fetch_span_for_tile` family: `ends`, `lengths`, `fragment-kmers`,
    - `midpoints`,
    - `ref-gc-bias`.
  - Do not use a reference that lacks the processed chromosome as the sentinel. A run that asks to
    process a chromosome absent from the reference is invalid and should fail.
  - Status: instrumentation/refactor-gated. The public command boundary currently has no honest
    way to distinguish "opened reference for setup" from "read sequence for a skipped tile" without
    making an invalid reference fixture. Add this when reference sequence loading can be injected,
    spied on, or tested at the tile helper boundary after pruning is moved before the read.

- Exact-output regression tests:
  - Add sparse-window fixtures where `fetch_span.start() > tile.fetch_start()`.
  - Use a reference sequence where shifting by the old tile origin would produce a different GC bin.
  - Use file-based GC correction with non-neutral weights so wrong GC changes output.
  - Assert parsed outputs exactly for:
    - `fcoverage`
    - `wps`
    - `midpoints`
    - `ends`
    - `lengths`
    - `fragment-kmers`
    - `gc-bias`
  - For `ref-gc-bias`, assert GC counts and ACGT totals match existing boundary-crossing behavior.
  - Add a command-level `gc-bias` regression where the owning tile has non-zero sequence origin and
    the saved observed count matrix has all mass at the hand-derived GC percentage.
  - Current sentinel fixture:
    - BAM chromosome length is 1,500 bp.
    - reference length is 1,022 bp, with `[0,900)` all A, `[900,961)` all C, and `[961,1022)` all A.
    - selected fragment is `[900,961)`.
    - GC package uses two bins, `[0,51)` weight 2.0 and `[51,101)` weight 7.0.
    - correct behavior reads only the narrowed span, converts `[900,961)` to prefix-local
      coordinates, computes 100% GC, and produces weight 7.0.
    - full-tile reads or old-origin indexing must either fail or produce the wrong output.

- Boundary and halo tests:
  - Include at least one chromosome-end case where the narrowed fetch span is clamped by chromosome
    end but still covers every eligible fragment.
  - Include at least one fragment whose aligned span overlaps the candidate window but whose
    contracted GC window is outside the loaded prefix span, to prove `neutralize_invalid_gc` and
    skip behavior are still intentional.
  - Status for contracted-GC-window outside-prefix command coverage: refactor-gated. Current
    commands either fetch a fragment only after window-derived pruning or use full-tile prefix
    spans, so this exact command-level state is created by the G-006 prefix narrowing refactor. The
    current `gc-bias::get_fragment_gc_returns_none_when_fragment_is_outside_loaded_sequence` helper
    test pins the intended `None` behavior until command-level coverage can be added.
  - Include blacklist masking with a narrowed reference origin, because masking must subtract
    `coordinates.reference_interval.start()` before mutating local sequence bytes.
  - Add a fragment-kmers regression that can distinguish the GC-correction reference slice from
    the k-mer extraction reference slice. The current late-tile sentinel proves the GC weight and
    k-mer output together, but it does not by itself prove those two loaded sequences cannot be
    accidentally conflated in a later rewrite.

## High-Risk Points

- Do not let the byte slice imply coordinates. It does not. Always pass a `GCPrefixCoordinateMap` when the slice came from a reference range.
- Do not use `tile.fetch_start()` as a prefix origin after narrowing. The prefix coordinate map owns the origin.
- Do not use `tile.fetch_start()` / `tile.fetch_end()` as GC-prefix coverage bounds after narrowing. The narrowed `fetch_span` is the loaded-prefix coverage.
- Do not name methods `gc_count_absolute`; the count is not absolute. The interval argument is in reference coordinates.
- Do not remove local prefix APIs. Some code intentionally uses prefix-local windows.
- Keep blacklist masking aligned with `coordinates.reference_interval.start()`.
- Keep `fragment-kmers` GC-correction sequence separate from its k-mer extraction sequence.
- Do not change window ids, group ids, sparse row ids, reducer keys, or output ordering.

## Verification Commands

Per repo rules, implementation should run:

```bash
cargo check
cargo check --tests
```

Do not run tests unless explicitly requested.
