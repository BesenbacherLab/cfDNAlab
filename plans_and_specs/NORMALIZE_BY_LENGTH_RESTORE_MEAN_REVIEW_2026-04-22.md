## `--normalize-by-length` restore-mean review of commit 937c5325

Date: 2026-04-22

Scope: commit `937c53251ee9da9128b02f9b3c577e60100b7aec` ("Adds restore-mean to
normalization mode"). This commit is the
current `ends` branch tip, so every change below is still present in the tree.
Reviewed against
[NORMALIZE_BY_LENGTH_RESTORE_MEAN_PLAN_2026-04-22.md](./NORMALIZE_BY_LENGTH_RESTORE_MEAN_PLAN_2026-04-22.md).

Legend: **[bug]** = likely incorrect behavior, **[design]** = correctness OK
but questionable architecture or usability, **[nit]** = cosmetic.

---

### What the commit matches in the plan

For context before the criticism. These parts look right:

- `LengthNormalizationMode { Off, UnitMass, RestoreMean }` as a separate enum
  in [config.rs:19-24](../src/commands/fcoverage/config.rs#L19-L24).
- `uses_length_normalization()` and `restores_mean_after_length_normalization()`
  helpers, replacing the old boolean checks cleanly.
- Bare `--normalize-by-length` maps to `unit-mass` via `default_missing_value = "unit-mass"`.
- Aligned `--by-size` fast path is disabled under restore-mean
  ([fcoverage.rs:1215](../src/commands/fcoverage/fcoverage.rs#L1215)) and falls
  back into the shared reducer path, as the plan requires.
- `scale_reduced_row` scales `coverage_sum` by `k` and `coverage_sum_of_squares`
  by `k^2`, which is mathematically correct for restoring a linear per-base
  scale.
- `mean_normalization_length` is reported in the run summary via
  `extra_statistics` for any length-normalization mode, not only restore-mean.
- `fragment-count-weights` still goes through `LengthNormalizationMode::UnitMass`
  with a specific test asserting this
  ([coverage_weights_tests.rs:394-416](../src/commands/coverage_weights/coverage_weights_tests.rs#L394-L416)).
- Warning path for `counted_fragments == 0` exists and keeps the output in
  unit-mass scale.
- Positional outputs get a real scaled-merge path
  (`merge_scaled_positional_tiles`) instead of silently concatenating raw
  per-tile files.

---

### [bug] 1. `LengthNormalizationMode::Off` is exposed to the CLI

In [config.rs:17-24](../src/commands/fcoverage/config.rs#L17-L24), the enum
derives `clap::ValueEnum`, which (by default) publishes every variant as a
valid value. That means users can pass
`--normalize-by-length=off` and get a silent no-op. The plan's user-facing
shape was: omit the flag entirely, bare flag = unit-mass, explicit
`--normalize-by-length=restore-mean`. `off` was only meant to exist internally.

**Fix:** mark the variant with `#[value(skip)]` (or `#[clap(skip)]`, depending
on the version) so it stays the internal default but is not selectable from
the CLI, OR rename to something users can't reach via kebab-case.

---

### [design] 2. Two coexisting setters (`_bool` and `_mode`)

[config.rs:342-355](../src/commands/fcoverage/config.rs#L342-L355) kept both:

- `set_normalize_by_length_bool(bool)` — legacy shim
- `set_normalize_by_length_mode(mode)` — the real setter

The `_bool` form is used only by tests
([fcoverage_tests.rs](../src/commands/fcoverage/fcoverage_tests.rs) and
[test_fcoverage_command.rs](../tests/test_fcoverage_command.rs)) and not by any
production call site, but it still covers ~12 test call sites.

`coverage_weights.rs:234-238` takes a `bool` parameter and maps it to
`UnitMass`/`Off` via an inline ternary — the shim would have fit here but wasn't
used.

**Fix:** migrate the test call sites to
`set_normalize_by_length_mode(LengthNormalizationMode::UnitMass)` and delete
the `_bool` shim. Either choose one pattern or remove the shim entirely.

---

### [bug] 3. Counter field names leak "restore_mean" into the UnitMass path

In [counters.rs:100-101](../src/commands/counters.rs#L100-L101):

```rust
restore_mean_tile_owned_fragments: u64,
restore_mean_tile_owned_normalization_length_sum: u64,
```

These fields are incremented whenever
`normalization_length_for_fragment` returns `Some(_)`, i.e. whenever
`uses_length_normalization()` is true. That includes plain `UnitMass`. The
fields are therefore populated in modes that do not restore anything.

The downstream computation in
[fcoverage.rs:576-585](../src/commands/fcoverage/fcoverage.rs#L576-L585) also
computes `mean_normalization_length` for both modes.

**Fix:** rename to `tile_owned_normalization_fragments` and
`tile_owned_normalization_length_sum`. This is load-bearing: the counter field
names are how someone reading the run log or downstream counter code decides
what the stat means, and the current names suggest "only filled in
restore-mean".

---

### [bug] 4. Silent divergence from the plan on the denominator accumulator

The plan says (section 3, 3b, and the test matrix):

> the sum must only increase when `was_counted` is true, matching the existing
> `counted_fragments` counter semantics

Codex implemented a *different* rule in
[fcoverage.rs:1410-1416](../src/commands/fcoverage/fcoverage.rs#L1410-L1416):

```rust
fn fragment_is_owned_by_tile_for_restore_mean_stats(fragment, tile) -> bool {
    let fragment_start = fragment.start();
    fragment_start >= tile.core_start() && fragment_start < tile.core_end()
}
```

This gates the accumulator on *tile ownership via fragment start*, not on
`was_counted`.

**Why I suspect codex did this:** the existing `counted_fragments` counter
double-counts fragments visible in multiple tile fetch halos (each tile that
overlaps the fragment adds its contribution to the core and sets `was_counted
= true`). If the denominator accumulator matched `counted_fragments`, it would
*also* double-count, inflating `sum_counted_length` relative to the true
fragment count.

**Why this is bad as written:**
- It silently diverges from the plan with no comment explaining the tradeoff.
- It decouples the denominator sample size from `counted_fragments`, so the
  two numbers in the run summary can disagree. There is no logging of the
  disagreement.
- It introduces an invariant ("fragment.start ∈ exactly one core") that holds
  for tile layouts we currently use but is nowhere asserted. A future tile
  layout (e.g. overlapping cores, fragments whose start sits in an unfetched
  gap between chromosomes) can silently miscount.
- If a fragment's start is *before the first tile's core* (e.g. chromosome
  edge case, read-without-pair situations), it contributes to coverage but
  never to the denominator. Silently excluded.

**Also note:** the comment at
[fcoverage.rs:1403-1408](../src/commands/fcoverage/fcoverage.rs#L1403-L1408)
warns the reader not to reuse this helper for generic fragment counting, but
does not explain *why it exists at all* — i.e. doesn't state the
multi-tile-halo double-counting concern that motivated it.

**Fix options, pick one and document it:**
1. Match the plan exactly — use `was_counted`, accept that the denominator may
   double-count, and document that the resulting mean is biased when halos
   cover fragments across tiles. (Simplest, matches counted_fragments.)
2. Keep the tile-ownership rule, but also change `counted_fragments` to the
   same rule so both stats stay consistent. (Bigger change, unbreaks
   counted_fragments if it's currently wrong.)
3. Keep the tile-ownership rule, rename the comment/doc to explain *why*, and
   add a debug assertion that `restore_mean_tile_owned_fragments <=
   counted_fragments` always holds.

---

### [design] 5. `tile_output_decimals = decimals_to_use.max(12)` is a hack

[fcoverage.rs:384-388](../src/commands/fcoverage/fcoverage.rs#L384-L388):

```rust
let tile_output_decimals = if opt.restores_mean_after_length_normalization() {
    decimals_to_use.max(12)
} else {
    decimals_to_use
};
```

This magic-constant bump exists because positional tile rows are later re-read
by `merge_scaled_positional_tiles`, multiplied, and rounded once more. If tile
rows were rounded to 0 or 2 decimals first, the re-scale would amplify
roundoff.

Problems:
- `12` is a magic constant with no explanation of why 12 (and not, say, 9 or
  15). For f64 you generally need ~15-17 significant digits to round-trip.
- It affects every tile-level write, even paths that use `writeln!("{}", f64)`
  for partials already. Those are already at full precision; the bump only
  truly matters for positional bedgraph/tsv tile files consumed by the scaled
  merger.
- The whole round-to-text-then-reparse round trip is avoidable. A raw f64
  intermediate (e.g. `ryu`, `{:?}`, or a binary sidecar) would sidestep both
  the precision choice and the CPU cost of re-parsing floats at merge time.

**Fix:** write positional tile values in a lossless text form (e.g. `ryu` or
`{:e}` with enough digits) and round only once at final merge. Delete
`tile_output_decimals`.

---

### [bug] 6. Zero-fragment warning condition doesn't match the plan intent

[fcoverage.rs:591-596](../src/commands/fcoverage/fcoverage.rs#L591-L596):

```rust
if opt.restores_mean_after_length_normalization() && mean_normalization_length.is_none() {
    warn!(...);
}
```

`mean_normalization_length.is_none()` evaluates to true iff
`restore_mean_tile_owned_fragments == 0`. Because of issue #4, that counter is
not the same as `counted_fragments`. So the warning can fail to fire when the
run counted fragments but none of them happened to be tile-owned by the
restore-mean rule, or fire erroneously in the mirror case.

Plan wording:

> If `counted_fragments == 0`, the mean is undefined.

Tie the warning (and the no-op behavior) to `counted_fragments == 0` after
aligning the accumulator semantics (issue #4).

---

### [nit] 7. Output filename suffix `length_normalized_with_restored_mean` is long

[fcoverage.rs:336](../src/commands/fcoverage/fcoverage.rs#L336) produces
filenames like:

```
<prefix>.fcoverage.length_normalized_with_restored_mean.per_position.bedgraph.zst
```

Consider a shorter, still-unambiguous suffix: `length_normalized.restored_mean`
or `length_normalized_restored`. Purely cosmetic but it will appear in every
output filename from now on.

---

### [design] 8. `scale_reduced_row` scales `coverage_sum_of_squares` even when unused

[writers.rs:346-359](../src/commands/fcoverage/writers.rs#L346-L359) always
scales both fields. The non-summary path
(`write_scaled_non_summary_reduced_row`) never reads
`coverage_sum_of_squares`, so the `powi(2)` multiply is dead work.

Harmless but noisy. Either:

- split into two helpers (sum-only / sum-and-SS), or
- leave a comment acknowledging the wasted multiply is intentional to keep one
  shared scale function.

---

### [design] 9. Grouped BED fold scales per segment, not per group

[writers.rs:680-735](../src/commands/fcoverage/writers.rs#L680-L735) scales
each reduced segment before folding into the group accumulator. This is
correct by linearity:

- `k*(A + B) == k*A + k*B` for `coverage_sum`
- `k^2*(SSA + SSB) == k^2*SSA + k^2*SSB` for `coverage_sum_of_squares`

So results are equivalent to scaling once at the end. But it multiplies `k` /
`k^2` N times per group instead of once. In tight reducer loops this may show
up in perf. Fine for correctness, possibly worth flipping to scale-at-write
later if you care about the cost.

---

### [design] 10. `sorted_tile_files_for_chromosome` refactor was bundled in

The commit extracts
[tiling.rs:21-46](../src/commands/fcoverage/tiling.rs#L21-L46) which
deduplicates the directory-scan/sort logic in three places. Good refactor, but
it's not a restore-mean change and should have been a separate commit. Easier
to audit, easier to revert restore-mean without reverting the good refactor.
Minor.

---

### [bug] 11. `merge_scaled_positional_tiles` preserves within-tile line order, not genomic order

[tiling.rs:106-217](../src/commands/fcoverage/tiling.rs#L106-L217) reads each
per-tile file in order and writes its lines in the order encountered. This
mirrors the non-scaled merger, so ordering behavior is consistent with the
existing path. But note the test at
[tiling_tests.rs:284-296](../src/commands/fcoverage/tiling_tests.rs#L284-L296)
(`merge_scaled_positional_tiles_orders_tiles_and_scales_values`) writes tile 0
as `["chr1\t0\t5\t...", "", "chr1\t10\t15\t..."]` — *out of internal order* —
and asserts the output preserves that internal order. That is probably fine
because real per-tile positional writers never emit out-of-order rows, but
the test is documenting "tile-internal order is preserved verbatim", not
"output is genomically sorted".

If a future writer ever emits positional rows in anything other than ascending
genomic order within a tile, the merged output silently becomes unsorted.
Worth adding either an invariant assertion or a test that validates genomic
ordering across tiles.

---

### [nit] 12. No sidecar file for `mean_normalization_length`

The plan marks this as optional. Not implemented. Fine to defer, but worth
tracking: if users ever need to reproduce a restore-mean rescaling offline (to
convert a unit-mass output into a restore-mean-equivalent output), the number
would need to be scraped from logs.

---

### [design] 13. `set_normalize_by_length_mode` takes the full enum including `Off`

Minor: the public API lets callers set `Off` explicitly. Combined with issue
#1, that's two paths for "no normalization". Not a correctness issue — just an
API surface observation.

---

## Summary: priority fixes

Ordered by impact / effort:

1. **#4 (denominator accumulator semantics)** — decide which rule is correct,
   document it, and make sure the accumulator and `counted_fragments` tell a
   consistent story. This is the one most likely to produce wrong numbers in
   production use.
2. **#3 (counter naming)** — low effort, high clarity.
3. **#1 (`Off` in CLI)** — one-line fix.
4. **#5 (tile_output_decimals hack)** — worth cleaning up before more paths
   accumulate magic-number precision bumps.
5. **#6 (warning condition)** — falls out once #4 is resolved.
6. **#2 (`_bool` shim cleanup)** — easy test-only migration.
7. **#8, #9** — cosmetic perf / clarity.
8. **#7, #10, #11, #12, #13** — nits.

---

## Failing-test f32-roundoff analysis

Example failure: `expected 120, got 120.000002` (test
`restore_mean_grouped_summary_stats_on_unique_bases_writes_scaled_rows`).

Why this happens:

- Coverage arrays are `f32` internally
  ([coverage.rs:736](../src/shared/coverage.rs#L736): `coverage() -> Option<&[f32]>`).
- Unit-mass per-base weights are `1 / counted_length`. For the test fixture
  (fragment 1: length 40, fragment 2: length 80), those weights are `0.025`
  and `0.0125`, neither of which is exactly representable in `f32`.
  - `0.025_f32` ≈ `0.025000000372529`
  - `0.0125_f32` ≈ `0.012500000186264`
- Sum of 40 copies of `0.025_f32` ≈ `1.000000014901`; sum of 80 copies of
  `0.0125_f32` ≈ `1.000000014901`; combined unit-mass `coverage_sum` ≈
  `2.000000029802`.
- `mean_normalization_length = (40 + 80) / 2 = 60`.
- After scaling: `2.0000000298 * 60 ≈ 120.00000179`, which rounds to
  `120.000002` at 6 decimals.
- Test tolerance is `1e-9`, which cannot absorb the amplified f32 error. The
  error is small (`~1.5e-6`) but real.

**This is not a bug introduced by codex.** The same f32 error is present under
unit-mass — it just wasn't visible at the test tolerance because the scale is
smaller (`2.0` vs `120.0`). Restore-mean amplifies the same relative error by
the multiplier, making it visible.

### Options (no code changes yet)

1. **Loosen test tolerance for restore-mean.** Pragmatic, minimal. The absolute
   tolerance should scale with `mean_normalization_length`. A tolerance of
   e.g. `1e-4 * mean_normalization_length` (or equivalently `1e-4` relative)
   would comfortably cover the observed ~`1.5e-6` absolute error while still
   catching real regressions. This is the path I'd recommend for now; the
   current `1e-9` tolerance was always pretending f32 coverage didn't exist.

2. **Compensated / Kahan summation in the prefix sum builder.** Reduces the
   summation error but doesn't fix the per-element f32 quantization, so you'd
   still see ~one ULP of f32 per element. Marginal benefit, real cost.

3. **Move the internal coverage array to f64.** Kills this class of error
   entirely. Doubles the memory cost of the coverage buffer. Given cfDNAlab's
   typical tile size, may or may not matter — worth measuring. Cleanest
   long-term fix.

4. **Apply restore-mean scaling during accumulation instead of at the end.**
   Means multiplying base_weight by `mean_normalization_length` before adding
   to the coverage array. This restores the signal to O(1) per base instead of
   O(1/length), but requires knowing `mean_normalization_length` during the
   counting pass, which we don't. The plan explicitly rejects this for exactly
   that reason. Not viable.

5. **Round tile-level and final outputs to fewer decimals.** Currently the
   failing test sets `decimals = 6`. Dropping to `decimals = 4` would hide this
   error but is a user-facing change to output precision, not a real fix.

6. **Choose test fixtures whose weights are f32-exact** (e.g. fragment lengths
   that are powers of two, so `1/length` has an exact f32 representation).
   Only fixes tests, not real behavior. Avoid as a general strategy.

7. **Track and subtract the known f32 bias analytically** in tests. Brittle and
   couples tests to the internals.

**Recommended path:**

- Short term: (1). Loosen the tolerance for restore-mean tests to `1e-4`
  relative, or `1e-4 * mean_normalization_length` absolute. Document the
  reason in a single shared test helper (`assert_close_for_restore_mean` or
  similar).
- Medium term: consider (3). f64 coverage arrays would also make other
  existing-but-hidden f32 precision issues go away. Worth evaluating the
  memory impact against current tile sizes.
- Avoid (2), (4), (5), (6), (7).
