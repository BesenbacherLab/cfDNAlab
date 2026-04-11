# `ends` Command — Deep Code Review

Scope: all files under `src/commands/ends/`, the shared `ends_fragment.rs`, and related
infrastructure called exclusively from `ends`. All claims are cross-referenced with the
actual source lines.

---

## Part 1 — Confirmed Bugs

### M11 — Statistics and progress go to `stdout`, not `stderr` (convention note)

**File:** `ends.rs:87, 99, 119, 163, 211, 227, 338-374`

All `println!()` calls (start banners, tile progress banner, statistics block) go to
stdout. Since `ends` writes all results to files and emits nothing useful on stdout,
there is no practical conflict today. The concern is purely conventional: stderr is the
standard channel for diagnostics so that scripts capturing stdout get a clean stream.
Not a bug in the current design — but worth tracking if stdout data output is ever
added.

---

### B4 — `keep_temp = false` is hardcoded dead code with an unreachable branch

**File:** `ends.rs:310-321`

```rust
let keep_temp = false;
if !keep_temp {
    if let Err(e) = std::fs::remove_dir_all(&temp_dir) { … }
} else {
    eprintln!("kept temp tiles in {}", temp_dir.display()); // ← unreachable
}
```

The `else` branch can never execute. This is a debugging leftover. Either expose
`keep_temp` as a config option (or env-var) or delete the `else` branch entirely.

---

## Part 2 — Logic / Semantic Issues

### L1 — `fragment_assignment_length` computed unconditionally but only used by `Midpoint`

**File:** `ends.rs:587`

```rust
let fragment_assignment_length = fragment.assignment_len();
```

This is computed for every fragment but only consumed inside the `Midpoint` arm of
the `query_interval` match (line 613). All other arms ignore it. Minor, but it is
noise that suggests the intent was to use it more broadly.

---

### L2 — `min_overlap_fraction` formula does not account for `ClipStrategy::RawShiftedBoundary` expansion beyond `max_fragment_length`

**File:** `ends.rs:490-500`

```rust
1. / (2. * opt.fragment_lengths.max_fragment_length as f64 + 1.0)
// 2x to allow for raw-clipping-mode expansion
```

For `ClipStrategy::RawShiftedBoundary`, each end's assignment boundary extends outward by
`soft_clip_bp`. If `max_soft_clips` is not set there is no bound on this expansion. The
`2×` multiplier assumes the expansion is at most `max_fragment_length`, which is the
aligned length limit, not the assignment interval limit. For extremely soft-clipped reads
(unusual but valid), the formula underestimates the required denominator and could miss
windows that do contain the assignment boundary.

---

### L3 — `counted_fragments` counts fragments that pass GC+window filtering but may still have no motifs emitted

**File:** `ends.rs:661`

```rust
counter.base.counted_fragments += 1;   // incremented here
```

This is incremented after the GC check (which may have `continue`d) but before the
actual `count_fragment_in_window` loop. A fragment that reaches this line will always
produce at least one `count_fragment_in_window` call, but both ends can individually
produce `None` from `maybe_encode_end_motif_key` (blacklisted reference bases, indels,
N bases). Such a fragment is counted in the statistics as "counted" but contributes
zero actual counts to any window. The label "Fragments counted one or more times" is
therefore slightly inaccurate.

---

### L4 — Double-building of `within_spec` / `outside_spec` per tile

**File:** `ends.rs:229-230` and `motifs.rs:124-125`

`build_optional_kmer_spec` is called once in `run()` (for the reduction phase) and
again inside `build_tile_motif_context` for every tile. The specs are cheap but the
duplication is confusing: a future change to one site will silently diverge from the
other.

---

### L5 — `TileMotifContext.ref_2bit` is `Option<PathBuf>` (clone of path); opens a new 2bit handle per fallback call

**File:** `motifs.rs:37-38, 612-625`

`fetch_reference_kmer_exact` opens the 2bit file fresh for every individual fallback
lookup. For boundary fragments with large `k_outside`, this means repeated file opens
per fragment. The path is cloned cheaply but the I/O is not. A per-tile 2bit handle
stored in `TileMotifContext` would avoid this.

---

### L6 — `encode_blacklist_validation_inside_code` returns `Some(sentinel_none)` for right end below `k`

**File:** `motifs.rs:462-465`

```rust
if (end.boundary_pos as u64) < k {
    return Ok(Some(spec.sentinel_none()));
}
```

The sentinel_none is returned, then checked by `motif_code_is_invalid` and causes the
end to be skipped. But the check in `encode_inside_code` (line 402) only acts on this
result if `motif_code_is_invalid` is true. The flow is:

```
encode_blacklist_validation_inside_code → Some(sentinel_none)
→ motif_code_is_invalid → true
→ return Ok(masked_reference_code)   // which IS the sentinel
```

Then the caller `maybe_encode_end_motif_key` sees the sentinel and returns `None`. This
is correct but convoluted. The two code paths for "position underflows" and "blacklisted
base" are merged into the same sentinel code and it is hard to trace which reason caused
the skip.

---

### M12 — `debug_assert_eq!` in `format_end_motif_label` is a sanity check, not a real guard

**File:** `counting.rs:192`

```rust
debug_assert_eq!(full_motif.len(), within_len + outside_len);
```

The invariant is structurally guaranteed: `decode_kmer` always returns exactly `k`
characters (Ns for sentinels included), so the concatenated `full_motif` length is
always `k_inside + k_outside`. The `debug_assert` would only fire if the codec was
changed to return variable-length strings, which is not a realistic failure mode. It
is fine as-is — just noting it is a sanity check rather than real protection.

---

### L9 — `build_all_end_motif_order` size guard uses `4^(k_inside + k_outside)` (ignores collapse reduction)

**File:** `output.rs:182-184`

```rust
let motif_count_upper = 4_u64.checked_pow(total_k as u32)…;
```

With `--collapse-complement`, the actual column count is ≈ `4^k / 2`. The guard is
conservative, which means it will reject some `(k, n_windows)` combinations that would
actually fit. Not a correctness issue, but users can hit the guard for parameters that
would work with collapse enabled. The guard could be tightened to report the actual
post-collapse estimate.

### L11 — `assignment_interval` in `fragment_with_two_ends` test helper is identical to `interval`

**File:** `motifs_tests.rs:66-67`

```rust
assignment_interval: Interval::new(left_boundary_pos, right_boundary_pos)
    .expect("valid assignment interval"),
```

This sets `assignment_interval == interval` for all test fragments. In
`ClipStrategy::RawShiftedBoundary`, these would differ. Since the tests use `KmerSource::Read`, the
assignment_interval is only used for midpoint/overlap logic inside `process_tile` (not
in `count_fragment_in_window` itself), so the tests are still valid — but the helper is
misleading for anyone writing future tests involving Raw clip mode or non-Endpoint
assignment.

---

### L12 — `EndReadInfo::from_record_with_gc_tag` passes `clip_strategy` to `motif_has_indels` but `ClipStrategy::Aligned` and `ClipStrategy::Skip` are treated identically

**File:** `ends_fragment.rs:377-379`

```rust
let aligned_bases_in_motif = match clip_strategy {
    ClipStrategy::Aligned | ClipStrategy::Skip => k_inside,
    ClipStrategy::RawAlignedBoundary | ClipStrategy::RawShiftedBoundary => {
        k_inside.saturating_sub(soft_clip_bp as usize)
    }
};
```

`Skip` and `Aligned` behave identically in the indel check. But `Skip` mode later
skips the end entirely if `soft_clip_bp > 0` (line 476). The indel check for a
soft-clipped `Skip` end is therefore wasted work. Not a bug, but an unnecessary
computation that the structure implies may need special-casing if `Skip` semantics
ever change.

---

### L13 — `collect_fragment_with_ends_from_single_read` calls `resolve_fragment_end` on the *same* `read` for both ends

**File:** `ends_fragment.rs:268-322`

For `--reads-are-fragments`, both the left and right ends are resolved from the same
`EndReadInfo` object. The `left_soft_clip_bp` / `right_soft_clip_bp` are defined per-
read and should correctly match the left and right ends. But `is_inwards_oriented` is
never called, so a read that is in a non-inward orientation would still produce a
fragment. For single-read mode, orientation is not meaningful, so this is probably
intentional — but it diverges from the paired-end flow in a way that is undocumented.

---

## Part 3 — Configuration & API Inconsistencies

### C1 — `EndsConfig` programmatic API setters have no cross-field validation

**File:** `config.rs:288-330`

`set_ref_2bit(None)` after `source_inside = KmerSource::Reference` leaves an
inconsistent state. `run()` detects this at tile time (line 443 in ends.rs), but the
error is deferred until the processing begins, not caught at configuration time.

---

### C3 — `WindowMotifAssigner::FromStr` and clap `ignore_case = true` are inconsistent for direct `FromStr` callers

**File:** `config_structs.rs:155-180`, `config_structs.rs:116-120`

The clap attribute `ignore_case = true` causes clap to lowercase the input before
calling `FromStr`. Direct `FromStr` calls (in tests, downstream crates, or
deserialization) are case-sensitive. A caller passing `"Endpoint"` via `FromStr` gets
an error, while the CLI accepts it silently. Same issue applies to
`IndelMotifFilterPolicy` (config.rs:148).

---

### C6 — Hand-rolled JSON in `write_end_settings_json` is brittle

**File:** `write.rs:96-155`

JSON is manually assembled with `writeln!` and hardcoded commas. The final field
(`max_fragment_length`, line 151) correctly has no trailing comma, but the structure
is fragile: adding a new field requires correctly placing a comma on the
previously-last field. The project already depends on `serde`, making `serde_json` an
obvious alternative. If that dependency is to be avoided, at minimum a helper that
builds a `Vec<(key, value)>` and emits commas only between entries would be safer.

---

### C8 — `EndsConfig` does not implement `Default`; the TODO in `config_structs.rs:6` is stale

**File:** `config_structs.rs:6`

```rust
// TODO these structs should use the format used by other cfDNAlab commands instead
```

This TODO has presumably been open since the command was started. Without `Default`,
test setup and downstream programmatic use must call `EndsConfig::new()` or construct
all fields manually. The inconsistency with other commands creates friction.

---

## Part 4 — Test Coverage Gaps & Soundness Issues

### T12 — `decode_end_motif_counts_drops_motifs_with_n` constructs an N-containing key directly

**File:** `counting.rs:388-404`

```rust
inside_code: within_spec.encode_kmer_bytes(b"AN"),
```

`encode_kmer_bytes(b"AN")` returns `sentinel_n` (not a code for 'A'+'N', since 'N'
maps to digit 4 which triggers sentinel). The test is therefore testing that
`sentinel_n`-encoded keys produce N-strings that are then dropped. This is the correct
path. However, the test comment "sentinel-N decodes to an N-containing motif" only
partially explains the mechanism — it would be clearer to say "encode_kmer_bytes with
any N base returns sentinel_n, which decodes to an all-N string".

---

## Part 5 — Performance Observations

### P1 — `KmerSpec` rebuilt per tile inside `build_tile_motif_context`

Already noted in L4. In a run with 1000 tiles, `build_kmer_specs` is called ~2000
times total (once in `run()` + once per tile). The cost is trivial but the duplication
is a maintenance hazard.

---

### P2 — `scaling_with_bin_idx` is rebuilt from scratch for every tile

**File:** `ends.rs:510-513`

```rust
let scaling_with_bin_idx: Vec<IndexedInterval<u64>> = scaling_chr
    .iter()
    .map(|(start, end, _)| IndexedInterval::new(*start, *end, 0_u64))
    .collect::<crate::Result<_>>()?;
```

For each of N tiles on the same chromosome, this allocates a new `Vec` wrapping the
same scaling intervals. This could be precomputed once per chromosome outside the tile
loop (before `par_iter`) and referenced via a slice.

---

### P3 — `fragment_assignment_length` unnecessary allocation in non-Midpoint modes

**File:** `ends.rs:587`

Noted in L1. The result is only used in one arm. Inlining it into the Midpoint arm
would make the dependency explicit.

---

### P4 — `TileResult` carries `chr: String` but chr is derivable from the tile

**File:** `tiling.rs:26`

`TileResult` stores a `String` for the chromosome name. This is used nowhere in the
reduction phase (only the `counts_path` and `counter` are read). The field is dead
weight in the struct but causes a heap allocation per tile.

---

## Part 6 — Documentation Issues

### D2 — `Endpoint` mode doc does not explain how it interacts with `CountOverlap` weight

**File:** `config_structs.rs:98-113`

The `--assign-by` documentation explains `count-overlap` as "count up the fraction of
fragment bases overlapping each window" but does not clarify that, in this mode, BOTH
ends of a fragment are counted in EVERY candidate window (not just the window the end
falls in). This is a significant semantic difference that users familiar with
`endpoint` mode may miss.

---

### D4 — `--k-inside = 0` and `--k-outside = 0` check is documented but the allowed combination (`k_inside=0 OR k_outside=0`) is not

**File:** `ends.rs:74-76`

The guard allows one of the two to be zero. Using `k_inside=0` means only outside
bases are counted; using `k_outside=0` means only inside bases. The motif label format
becomes `"GT_"` or `"_AC"`. This is not documented in the help text for either
`--k-inside` or `--k-outside`.

---

### D5 — `ClipStrategy::Skip` vs `--max-soft-clips 0` equivalence is not documented

Both cause soft-clipped ends to be skipped, and `Skip` uses `SkipEndDropAssignmentBoundary`
while `max_soft_clips = 0` also uses `SkipEndDropAssignmentBoundary`. For counting
purposes they are equivalent. This should be noted so users don't combine them
redundantly (or misunderstand that they're not subtly different).

---

### D6 — The `inside_bases` orientation convention for `KmerSource::Read` on reverse reads is undocumented

**File:** `ends_fragment.rs:560-573`

`extract_right_inside_bases` takes the last `k_inside` bases of `read.seq`. This relies
on the BAM convention that `seq` is stored in the forward-strand orientation for
reverse-strand reads (i.e., BAM stores the reverse-complement of the original read
sequence for reverse reads). This is a non-obvious assumption that is not documented in
any comment near the extraction function. A future reader might incorrectly assume `seq`
is always in read orientation and introduce a bug.

---

### D7 — Motif orientation convention (5'→3' from end inward) is documented in `config.rs` but not in `counting.rs` where the key encoding actually happens

**File:** `config.rs:177-179` (good doc), `counting.rs` (no such comment)

`decode_full_motif` implements the 5'→3' orientation by reverse-complementing right-end
motifs. The function comment (lines 198-215) explains the storage/decode sequence, but
does not state the final biological interpretation. A reader understanding the
intermediate steps could still misidentify what "outside" and "inside" mean
biologically after the RC.

---

## Part 7 — Minor / Style Issues

### M1 — `GCCorrector` is `clone()`d once per tile even though it is likely already wrapped in `Arc`

**File:** `ends.rs:199`

`gc_corrector.clone()` inside `par_iter`. If `GCCorrector` is not internally
`Arc`-wrapped, each tile gets its own full copy of the correction matrix. The comment
at line 118 says "Load GC correction matrix" suggesting it could be large. An explicit
`Arc<GCCorrector>` in the outer type would make clone costs explicit.

---

### M2 — Dead drop calls at end of `run()`

**File:** `ends.rs:213-219`

```rust
drop(chr_offsets_for_threads);
drop(tile_window_spans_for_threads);
drop(tile_window_spans);
drop(tiles);
drop(scaling_map);
drop(gc_corrector);
```

These explicit drops are unnecessary: Rust drops values when they go out of scope.
Explicit drops only have a semantic effect when you need to release something *before*
a subsequent blocking operation (e.g., releasing a mutex before waiting). Here they
precede the reduction phase but the comment just says "Release per-tile inputs before
merging outputs" — there is no locking involved. Remove or document the actual
motivation.

---

### M3 — `inspect_cigar_edges` would silently accumulate multiple consecutive soft-clip operations

**File:** `ends_fragment.rs:342-366`

SAM spec forbids two consecutive `S` operations, but the loop would add them both if
they appeared. The behavior is technically correct for invalid CIGAR (garbage-in
garbage-out), but a debug_assert or comment noting this assumption would help.

---

### M4 — `ClipStrategy` derives `Default` (Aligned) but `ClippingArgs` does not derive `Default`

**File:** `config_structs.rs:31-80`

`ClipStrategy` has `#[default]` on `Aligned`. `ClippingArgs` does not derive
`Default`, so it cannot participate in struct update syntax (`..Default::default()`).
This is inconsistent: if the inner enum has a meaningful default, the struct wrapping
it should too.

---

### M5 — `TileResult.chr` field is unused after construction

**File:** `tiling.rs:26`, `ends.rs:756-760`

```rust
TileResult {
    chr: tile.chr.clone(),  // stored but never read
    counts_path,
    counter,
}
```

`chr` is set from `tile.chr` but neither the counter accumulation loop nor the merge
loop ever reads `tile_result.chr`. It's dead data.

---

### M6 — `split` over `&[prefix, "counts"]` to build `counts_prefix` creates a string on every tile naming

**File:** `ends.rs:161`

```rust
let counts_prefix = &dot_join(&[prefix, "counts"]);
```

`dot_join` allocates a `String`. That string is then used as a template for every tile
file name (line 433-437). This is fine — it's a one-time allocation — but the
`counts_prefix` is a `&str` borrowed from a temporary. In Rust this compiles because
the temporary lives to the end of the function, but the notation (`&dot_join(...)`)
binding a reference to a temporary feels fragile and might confuse readers unfamiliar
with Rust's lifetime rules.

---

### M7 — `EndsConfig::new()` silently accepts `k_inside=0` and `k_outside=0` simultaneously

**File:** `config.rs:246-286`, `ends.rs:74-76`

The validation `k_inside == 0 && k_outside == 0` happens only in `run()`, not in
`new()`. Constructing a config with both zero is possible programmatically without
an error until `run()` is called.

---

### M8 — `format_end_motif_label` called `format_end_motif_label` but just splits and rejoins with underscore

**File:** `counting.rs:185-196`

The function name is accurate but its docstring says "The full motif string is expected
to be oriented already and ordered as `outside || inside`." This ordering assumption
should be an explicit pre-condition (assert or type-level guarantee), not just a
docstring note.

---

### M9 — `motifs.rs` `EndSide` enum is `pub(crate)` but could be local to the module

**File:** `motifs.rs:62`

`EndSide` is only used within `motifs.rs`. Making it `pub(crate)` leaks the abstraction
unnecessarily. The `pub(crate)` on `TileMotifContext` and `count_fragment_in_window` is
justified (they are called from `ends.rs`), but `EndSide` is never used outside this
file.

---

### M10 — `maybe_collapse_full_motif` takes `motif: String` by value (forces clone at call site)

**File:** `counting.rs:256`

```rust
pub fn maybe_collapse_full_motif(motif: String, collapse_complement: bool) -> String {
```

The caller `decode_full_motif` (line 154) passes an owned `String` so no clone is
needed there. But the call in `build_all_end_motif_order` (line 75) constructs a
temporary `String` just to pass in:

```rust
format!("{outside}{inside}")
```

This is fine as-is but the function signature forces the caller to relinquish ownership
even when the result is discarded, which is counterintuitive.

---

## Part 8 — Integration Test Deep Analysis

Scope: `tests/test_ends_command.rs` (56 tests) plus the shared `tests/fixtures/mod.rs`.
All line numbers are from the integration test file. The analysis covers soundness of
assertions, correctness of mental derivations in comments, assertion gaps, and structural
issues.

---

### 8.1 Test Fixture Limitations

#### IT-I1 — `paired_fragment` fixture produces identical left and right motif labels under `KmerSource::Read`

**File:** `fixtures/mod.rs:59,72`

`paired_fragment(start, fragment_len, read_len)` creates:
- Forward read: `seq = vec![b'A'; read_len]`
- Reverse read: `seq = vec![b'T'; read_len]`, `is_reverse = true`

For `k_inside=1, k_outside=0, KmerSource::Read`:
- Left inside = first 1 base of forward seq = `'A'` → label `_A`
- Right inside = last 1 base of reverse seq = `'T'`, decoded as `RC('T') = 'A'` → label `_A`

Both endpoints produce the **same motif label**. Any integration test using this
fixture with `KmerSource::Read` cannot distinguish left from right endpoint by motif
label. A regression that swapped left/right endpoints would pass every such test as
long as the counts sum correctly.

The two CLI tests do not pass `--source-inside`, so they use the Read default and both
endpoints give `_A`. The assertion `"Distinct counted end motifs: 2"` in the tile-skip
test (line 422) counts 2 endpoint *instances*, not 2 distinct motif *strings*. The stat
field name "Distinct" is misleading in this case.

---

#### IT-I2 — `simple_reference_twobit` is a strictly periodic `ACGT` repeat

**File:** `fixtures/mod.rs:180`

```rust
let chr1 = ("chr1".to_string(), "ACGT".repeat(64));
```

`chr1` is 256 bases of `ACGT` repeating. Every reference position is fully
deterministic: `base(i) = "ACGT"[i % 4]`. The tests rely on this to drive specific
expected values (e.g., `reference[10] = 'G'`, `reference[19] = 'T'`). This is good
for exact testing, but:

- Outside-k tests with `k_outside ≥ 4` will see repeating patterns. The fallback
  test (`outside_reference_lookup_falls_back_to_exact_reference_fetch`, line 994) uses
  `k_outside=5` and correctly derives `"CGTAC"` from positions 5–9 using the period-4
  sequence. If the reference were non-periodic, that derivation would differ. The test
  does not guard against accidentally correct results from the periodicity.

- A reference with `N` bases at specific positions, which would test the `sentinel_n`
  path in motif encoding, is only constructed in `skip_invalid_gc` (line 785) for GC
  purposes, not for motif encoding.

---

#### IT-I4 — `single_read_bam` produces reads with `flags: 0`

**File:** `test_ends_command.rs:97`

```rust
flags: 0,
```

Flags = 0 means: not paired, not properly paired, not reverse strand, not first in
pair, etc. This produces reads that real single-end sequencing would not emit — real
single reads have the `READ_UNMAPPED`/`READ_PAIRED` flags differently. The unpaired
tests use `reads_are_fragments = true`, which bypasses pairing logic, so the flags are
irrelevant in practice. But a test reader might wonder why single reads have flags = 0.

---

### 8.2 Assertions That Miss Regressions

#### IT-S3 — `midpoint_assigns_both_end_motifs_to_the_midpoint_window` comment is imprecise

**File:** `test_ends_command.rs:1084`

```
// fragment [10,20) has even midpoint 14 or 15, both inside [14,16).
```

Fragment start=10, aligned_length=10. The midpoint formula in `ends.rs` is:
`start + fragment_assignment_length / 2 = 10 + 10/2 = 10 + 5 = 15` (integer division).

The midpoint is **deterministically 15**, not "14 or 15". The comment acknowledges
ambiguity that does not exist in the code. If the formula ever changed to ceiling
division (midpoint = 15), the test would still pass since the window is `[14,16)`.
The comment should say "midpoint 15", not "14 or 15".

---

#### IT-S4 — Float equality assertions rely on IEEE 754 exact representability

**File:** `test_ends_command.rs:1153-1154, 1269-1271`

```rust
assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 0.5);  // count_overlap
assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);  // proportion
```

The `count_overlap` test uses window `[10,15)` and fragment `[10,20)`: overlap = 5/10
= 0.5. The value 0.5 = 1/2 is exactly representable in f64, so `5.0_f64 / 10.0 ==
0.5` is guaranteed true. If a future test used overlap 6/10 = 0.6 (not exactly
representable), the same `assert_eq!` pattern would be fragile. The suite should
document the IEEE 754 dependency or switch to approximate comparisons for non-dyadic
fractions. The GC-scaling combined test (`blacklist_gc_and_scaling_weights_combine_*`)
uses weights like `2.5` and `0.75` — both dyadic, so exact. Currently safe but
undocumented.

---

#### IT-S5 — `proportion_assignment` exact-boundary test is actually a good positive case

**File:** `test_ends_command.rs:1256`

```rust
assign_by: WindowMotifAssigner::Proportion(0.5),
```

Fragment `[10,20)`, window `[10,15)`, overlap = 5/10 = 0.5, threshold = 0.5. The
comparison `0.5 >= 0.5` is true, so the fragment is accepted. This is a correct
boundary test. The complementary rejection test (line 1278) uses window `[10,14)`,
overlap = 4/10 = 0.4, which is unambiguously below 0.5. The pair is sound. ✓

The edge case `proportion = 0.0` (from L10) has no test, so the behavior that every
fragment is accepted when threshold is 0.0 is untested.

---

#### IT-S6 — `cli_statistics_only_count_reads_from_tiles_with_relevant_windows` stat interpretation is ambiguous

**File:** `test_ends_command.rs:422`

```rust
assert!(stdout.contains("Distinct counted end motifs across those fragments: 2"));
```

Fragment B uses `KmerSource::Read` (no `--source-inside` flag). Forward seq = AAAA,
reverse seq = TTTT. Left inside = 'A' → `_A`. Right inside = RC('T') = 'A' → `_A`.
Both endpoints emit the same motif label `_A`, so there is 1 distinct motif *string*
but 2 endpoint counting *events*. The counter reports 2, which the test asserts.

The stat label "Distinct counted end motifs" appears to mean "endpoint instances
counted", not "distinct motif labels". A user reading the printed statistics might
interpret "distinct" as the number of unique motif labels and be confused when a
fragment with 2 identical endpoints reports "2 distinct end motifs". The test
inadvertently documents the confusing stat name.

---

#### IT-S7 — `raw_endpoint_assignment` comment derivation is correct but subtle

**File:** `test_ends_command.rs:1579-1581`

```
// The raw terminal bases are T on the left and A on the right, which both orient to "_T".
```

- Left raw base: `seq[0]` of `b"TTACGTAA"` = `'T'`. Left inside = `'T'` → label `_T`. ✓
- Right raw base: `seq[7]` = `'A'`. Right inside = `'A'`, decoded as `RC('A') = 'T'` → label `_T`. ✓

The comment correctly notes that both ends produce `_T` for different reasons (T on
the left is kept as-is; A on the right is reverse-complemented). Sound.

---

### 8.3 Mental Derivation Checks

#### IT-C1 — `outside_reference_lookup_falls_back_to_exact_reference_fetch` derivation is correct and rigorous

**File:** `test_ends_command.rs:996-1000`

```
// ...Asking for k_outside=5 on the left endpoint at 10
// needs reference bases [5,10), which crosses one base outside the preloaded tile slice and
// must therefore use the exact fallback path. On the ACGT-repeat reference, seq[5..10) is
// C G T A C, so the outside-only label is "CGTAC_".
```

Verification:
- tile_size=100 (default in `base_config`), max_frag_len=4
- Tile halo = max_frag_len = 4. Tile 0 starts at 0; local slice loaded with extra 4
  bp before pos 0... but this is `[10,11)` window so the tile covers [0, 100). The
  halo for k_outside prefetching extends the reference slice back by k_outside from the
  leftmost window boundary. But the *exact fallback* is triggered when k_outside exceeds
  max_frag_len. Here k_outside=5 > max_frag_len=4, so the tile slice covers [6, …) but
  position 5 is needed → fallback. Correct reasoning. ✓
- `reference[5..10]`: pos 5=C, 6=G, 7=T, 8=A, 9=C → `"CGTAC"`. ✓
- Outside-only label format: `"CGTAC_"` (outside bases then underscore for empty inside). ✓

This test has the most rigorous mental derivation in the suite.

---

#### IT-C2 — `blacklist_masking_still_skips_read_backed_within_motifs` derivation note

**File:** `test_ends_command.rs:430-432`

```
// - left read-backed motif = "_A"
// - right read-backed motif = reverse-complement("A") = "_T"
```

Read seq = `b"ACGA"`, so:
- Left inside = `seq[0] = 'A'` → `_A`. ✓
- Right inside = `seq[3] = 'A'`, decoded as `RC(b'A') = b'T'` → `_T`. ✓

The comment writes "reverse-complement("A") = '_T'", which applies RC to the *base*
`'A'` to get `'T'`, then forms the label `_T`. This is correct but could be read as
applying RC to the motif string `"A"` (which also gives `"T"` for a 1-mer). The
description is unambiguous for k=1 but would be wrong for k>1 if interpreted as
motif-level RC. Minor documentation risk for future test authors.

---

### 8.4 Missing Coverage

#### IT-M5 — Most tests use exactly one fragment

**File:** Most tests

Only three tests use multiple fragments: `cli_statistics_only_count_reads_from_tiles_with_relevant_windows`
(2 fragments), `blacklist_gc_and_scaling_weights_combine_*` (3 fragments), and
`cross_tile_fragment_is_counted_once_per_window_when_it_reaches_into_the_next_tile`
(1 fragment spanning 2 tiles). All window-assignment tests, all clip-mode tests, all
indel-filter tests, and all motif-label tests use exactly 1 fragment. Multi-fragment
window accumulation — where motif counts stack across fragments in the same window — is
essentially untested.

---

### 8.5 Structural Issues

#### IT-ST6 — All fixture fragments use the same start/length geometry

**File:** Most tests use `simple_paired_fragment_bam("...", 10, 10, 4)`

With rare exceptions (edge tests, cross-tile), every standard fragment is at position
10 with aligned length 10. Bugs that only manifest for fragments at position 0,
fragments whose k_outside window crosses the chromosome start, or fragments with
lengths that interact with tile-size modulo arithmetic would not be caught by the
standard geometry.

---

### 8.7 Test Suite Positive Observations

For completeness, the following tests are particularly thorough:

- **`blacklist_gc_and_scaling_weights_combine_to_the_exact_expected_endpoint_counts`**
  (line 552): Three fragments with a blacklisted region, separate GC correction
  weights per fragment, and per-region scaling factors, with exact per-motif
  numerical assertions. This is the highest-fidelity end-to-end numerical test.

- **`outside_reference_lookup_falls_back_to_exact_reference_fetch_when_the_motif_crosses_the_tile_slice`**
  (line 994): Complete mental derivation of the tile-slice halo boundary condition,
  verified against the exact reference sequence. The only test exercising the fallback
  reference lookup code path.

- **`cross_tile_fragment_is_counted_once_per_window_when_it_reaches_into_the_next_tile`**
  (line 1160): Validates the halo/deduplication invariant for a fragment spanning two
  tile boundaries. Critical behavioral guarantee with exact per-window assertions.

- **`raw_endpoint_assignment_uses_the_shifted_assignment_boundaries`** +
  **`aligned_endpoint_assignment_ignores_raw_shifted_boundary_positions`** +
  **`fragment_length_filters_use_the_aligned_fragment_length_even_in_raw_mode`**
  (lines 1575–1706): Three complementary tests that pin the Raw clip-strategy behavior
  from three different angles (assignment boundaries, window selection, length filter).

---

## Summary Table

| # | Category | Severity | File(s) |
|---|----------|----------|---------|
| M11 | Style | Low | ends.rs |
| B4 | Bug | Low | ends.rs |
| L1 | Logic | Low | ends.rs |
| L2 | Logic | Medium | ends.rs |
| L3 | Logic | Low | ends.rs |
| L4 | Logic | Low | ends.rs, motifs.rs |
| L5 | Logic | Low | motifs.rs |
| L6 | Logic | Low | motifs.rs |
| M12 | Style | Low | counting.rs |
| L9 | Logic | Low | output.rs |
| L11 | Logic | Low | motifs_tests.rs |
| L12 | Logic | Low | ends_fragment.rs |
| L13 | Logic | Low | ends_fragment.rs |
| C1 | Config | Low | config.rs |
| C3 | Config | Medium | config_structs.rs |
| C6 | Config | Low | write.rs |
| C8 | Config | Medium | config_structs.rs |
| T12 | Tests | Low | counting.rs |
| P1 | Perf | Low | motifs.rs |
| P2 | Perf | Low | ends.rs |
| P3 | Perf | Low | ends.rs |
| P4 | Perf | Low | tiling.rs |
| D2 | Docs | Medium | config_structs.rs |
| D4 | Docs | Medium | config.rs |
| D5 | Docs | Low | config_structs.rs |
| D6 | Docs | High | ends_fragment.rs |
| D7 | Docs | Medium | counting.rs |
| M1 | Style | Low | ends.rs |
| M2 | Style | Low | ends.rs |
| M3 | Style | Low | ends_fragment.rs |
| M4 | Style | Low | config_structs.rs |
| M5 | Style | Low | tiling.rs, ends.rs |
| M6 | Style | Low | ends.rs |
| M7 | Style | Low | config.rs |
| M8 | Style | Low | counting.rs |
| M9 | Style | Low | motifs.rs |
| M10 | Style | Low | counting.rs |
| IT-I1 | Tests/Infra | Medium | fixtures/mod.rs, test_ends_command.rs |
| IT-I2 | Tests/Infra | Low | fixtures/mod.rs |
| IT-I4 | Tests/Infra | Low | fixtures/mod.rs |
| IT-S3 | Tests | Low | test_ends_command.rs:1084 |
| IT-S4 | Tests | Low | test_ends_command.rs:1153 |
| IT-S6 | Tests | Low | test_ends_command.rs:422 |
| IT-M5 | Tests | Medium | (no test) |
| IT-ST6 | Tests | Low | test_ends_command.rs |
