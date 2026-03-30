# `ends` Command — Deep Code Review

Scope: all files under `src/commands/ends/`, the shared `ends_fragment.rs`, and related
infrastructure called exclusively from `ends`. All claims are cross-referenced with the
actual source lines.

---

## Part 1 — Confirmed Bugs

### B1 — `EndsConfig::new()` defaults `ClipStrategy::Raw`; CLI defaults `"aligned"`

**File:** `config.rs:259-262` vs `config_structs.rs:71`

```rust
// config.rs:259-262  (programmatic default)
clip: ClippingArgs {
    clip_strategy: ClipStrategy::Raw,
    max_soft_clips: None,
},

// config_structs.rs:71 (CLI default)
clap(long, value_enum, default_value = "aligned", …)
```

`ClipStrategy` derives `Default` and marks `Aligned` as the default
(`config_structs.rs:35`). The CLI correctly uses `aligned`. Every test or downstream
crate that calls `EndsConfig::new()` silently gets `Raw` instead. Any integration test
that builds config programmatically is testing a different code path than the user-facing
CLI.

---

### B2 — `EndsCounters::gc_out_of_range_tags` is declared but never incremented

**File:** `src/commands/counters.rs:127-131`, `ends.rs:645-659`

```rust
// counters.rs:127-131
counter_struct!(EndsCounters;
    blacklisted_fragments: u64,
    gc_failed_fragments: u64,
    gc_out_of_range_tags: u64   // ← never touched
);
```

In `ends.rs`, the GC match block increments only `gc_failed_fragments`.
`gc_out_of_range_tags` is never written. The field is dead weight, but because it is
`AddAssign`-accumulated across tiles, it will always read zero in the statistics. No
user-visible line in the statistics prints it, so it produces no observable artefact
today — but if someone ever adds a statistics line for it, it will silently report 0.

---

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

### D8 — Statistics only reflect reads in tiles that have relevant windows

**File:** `ends.rs:460-471`

When `fetch_span_for_tile` returns `None` (tile has no relevant windows), the function
returns before opening the BAM reader. This is correct behavior — scanning reads that
can't contribute to any output would be wasteful. However, the statistics block printed
at the end ("Total reads", "Initially accepted reads", etc.) only counts reads from
processed tiles. With a narrow BED window set, a significant fraction of reads is
never counted. The statistics header or a companion note should make this scope
explicit so users don't interpret the numbers as genome-wide totals.

---

### B6 — Scaling: `.context("no overlapping scaling bins found")?` panics on valid inputs

**File:** `ends.rs:675`

```rust
.context("no overlapping scaling bins found")?; // Should always find >= 1 bin
```

If the user provides a scaling file that does not cover every base of every chromosome
(gaps, chromosomes missing, etc.), a fragment in an uncovered region produces an opaque
`anyhow` error mid-run. The comment acknowledges the assumption but the error message
gives no guidance. The check should explain that scaling coverage is incomplete.

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

### L2 — `min_overlap_fraction` formula does not account for `ClipStrategy::Raw` expansion beyond `max_fragment_length`

**File:** `ends.rs:490-500`

```rust
1. / (2. * opt.fragment_lengths.max_fragment_length as f64 + 1.0)
// 2x to allow for raw-clipping-mode expansion
```

For `ClipStrategy::Raw`, each end's assignment boundary extends outward by
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

### L6 — `encode_blacklist_validation_within_code` returns `Some(sentinel_none)` for right end below `k`

**File:** `motifs.rs:462-465`

```rust
if (end.boundary_pos as u64) < k {
    return Ok(Some(spec.sentinel_none()));
}
```

The sentinel_none is returned, then checked by `motif_code_is_invalid` and causes the
end to be skipped. But the check in `encode_within_code` (line 402) only acts on this
result if `motif_code_is_invalid` is true. The flow is:

```
encode_blacklist_validation_within_code → Some(sentinel_none)
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
always `k_within + k_outside`. The `debug_assert` would only fire if the codec was
changed to return variable-length strings, which is not a realistic failure mode. It
is fine as-is — just noting it is a sanity check rather than real protection.

---

### L8 — Sparse path: `stack_end_motif_counts` silently drops motifs absent from `motif_columns`

**File:** `write.rs:187-190`

```rust
if let Some(&col) = motif_columns.get(motif) {
    mat[(row, col)] = count;
}
```

In dense `--all-motifs` mode, `motif_columns` is the full universe, so this branch is
dead. But any future code path that calls this function with a partial universe would
silently lose data. An assertion `assert!(motif_columns.contains_key(motif))` or an
explicit error would be safer.

---

### L9 — `build_all_end_motif_order` size guard uses `4^(k_within + k_outside)` (ignores collapse reduction)

**File:** `output.rs:182-184`

```rust
let motif_count_upper = 4_u64.checked_pow(total_k as u32)…;
```

With `--collapse-complement`, the actual column count is ≈ `4^k / 2`. The guard is
conservative, which means it will reject some `(k, n_windows)` combinations that would
actually fit. Not a correctness issue, but users can hit the guard for parameters that
would work with collapse enabled. The guard could be tightened to report the actual
post-collapse estimate.

---

### L10 — `proportion=0.0` is accepted as a valid window assigner

**File:** `config_structs.rs:172`

```rust
if !(0.0..=1.0).contains(&thr) { Err(…) }
```

A proportion of 0.0 means every window overlaps every fragment (since any overlap ≥ 0),
which would be pathological. This should probably be rejected or at least produce a
warning. There is no documentation indicating this is an intentional use case.

---

### L11 — `assignment_interval` in `fragment_with_two_ends` test helper is identical to `interval`

**File:** `motifs_tests.rs:66-67`

```rust
assignment_interval: Interval::new(left_boundary_pos, right_boundary_pos)
    .expect("valid assignment interval"),
```

This sets `assignment_interval == interval` for all test fragments. In
`ClipStrategy::Raw`, these would differ. Since the tests use `KmerSource::Read`, the
assignment_interval is only used for midpoint/overlap logic inside `process_tile` (not
in `count_fragment_in_window` itself), so the tests are still valid — but the helper is
misleading for anyone writing future tests involving Raw clip mode or non-Endpoint
assignment.

---

### L12 — `EndReadInfo::from_record_with_gc_tag` passes `clip_strategy` to `motif_has_indels` but `ClipStrategy::Aligned` and `ClipStrategy::Drop` are treated identically

**File:** `ends_fragment.rs:377-379`

```rust
let aligned_bases_in_motif = match clip_strategy {
    ClipStrategy::Aligned | ClipStrategy::Drop => k_within,
    ClipStrategy::Raw => k_within.saturating_sub(soft_clip_bp as usize),
};
```

`Drop` and `Aligned` behave identically in the indel check. But `Drop` mode later
skips the end entirely if `soft_clip_bp > 0` (line 476). The indel check for a
soft-clipped `Drop` end is therefore wasted work. Not a bug, but an unnecessary
computation that the structure implies may need special-casing if `Drop` semantics
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

### L14 — The `window_assigner_name` serialises `Proportion(f)` as `"proportion=<float>"` with no precision control

**File:** `write.rs:272`

```rust
WindowMotifAssigner::Proportion(value) => format!("proportion={value}"),
```

A user who passes `proportion=0.2` may get `"proportion=0.2"` back — or they may get
`"proportion=0.19999999999999998"` depending on how Rust formats the f64. The JSON
sidecar would then be misleading and could not be round-tripped through `FromStr`. The
precision should be fixed (e.g., `format!("proportion={value:.10}")` or the original
string preserved).

---

### L15 — `make_canonical` for odd-length motifs does not produce the same canonical as even-length (convention mismatch depends on parity of `k_within + k_outside`)

**File:** `src/shared/base.rs:68-78`

```rust
if mid == b'G' || mid == b'T' {
    return rev_complement(&kmer);
}
return kmer;
```

For even-length combined motifs (`k_within + k_outside` is even), `make_canonical`
returns the lexicographically smaller of the motif and its RC — the standard IUPAC
convention. For odd-length combined motifs the approach is a middle-base heuristic: if
the middle base is G or T, reverse-complement, otherwise keep.

These two rules can disagree. Concrete example with `k_within=2, k_outside=1`
(total k=3): motif `"GTA"` — RC is `"TAC"`. Lex min is `"GTA"` (G < T), but the
middle base is `T`, so the heuristic returns `"TAC"`. A user expecting lexicographic
canonicalization for all motifs would get the wrong representative.

This affects any run where `k_within + k_outside` is odd. Whether the middle-base
heuristic is the intended convention for cfDNA odd-length motifs should be documented
explicitly. If lexicographic canonicalization is desired uniformly, the even-length
branch should be applied for all lengths.

---

## Part 3 — Configuration & API Inconsistencies

### C1 — `EndsConfig` programmatic API setters have no cross-field validation

**File:** `config.rs:288-330`

`set_ref_2bit(None)` after `source_within = KmerSource::Reference` leaves an
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

### C4 — `fragment_length_basis` is hardcoded to `"aligned"` in the settings JSON

**File:** `write.rs:140`

```rust
writeln!(settings_writer, "  \"fragment_length_basis\": \"aligned\",")
```

This is a static string, not derived from `opt`. If the length basis is ever made
configurable (e.g., assignment-interval length), this will silently remain "aligned" in
all output files.

---

### C5 — Settings JSON omits many parameters that affect results

**File:** `write.rs:88-159`

The sidecar records: k_within, k_outside, source_within, clip_strategy, max_soft_clips,
indel_filter, window_assignment, collapse_complement, reads_are_fragments,
fragment_length_basis, min_fragment_length, max_fragment_length.

It does **not** record:
- `blacklist` files / `blacklist_strategy` / `blacklist_min_size`
- `scale_genome.scaling_factors`
- `gc.gc_file` presence / `gc.gc_tag`
- `min_mapq`
- `require_proper_pair`
- `all_motifs`
- `windows` spec (BED file path, window size, or global)

Reproducibility requires the full parameter set. A downstream user loading the `.json`
cannot reconstruct the exact run.

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

### C7 — `bins.bed` 4th column is `overlap_perc`, not a standard BED `name` column

**File:** `ends.rs:332`

```rust
writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
```

Standard BED4 uses the 4th column as `name`. Many downstream tools assume this. Using
it for a float percentage is non-standard and should be documented or the column
should be named something unambiguous (e.g., by adding a header, or using a 5-column
format with an explicit empty name).

---

### C8 — `EndsConfig` does not implement `Default`; the TODO in `config_structs.rs:6` is stale

**File:** `config_structs.rs:6`

```rust
// TODO these structs should use the format used by other cfDNAlab commands instead
```

This TODO has presumably been open since the command was started. Without `Default`,
test setup and downstream programmatic use must call `EndsConfig::new()` (which has the
`ClipStrategy::Raw` bug, B1) or construct all fields manually. The inconsistency with
other commands creates friction.

---

## Part 4 — Test Coverage Gaps & Soundness Issues

### T1 — `count_fragment_in_window_any_counts_both_ends_in_same_window` does not verify the actual keys or weights

**File:** `motifs_tests.rs:159-182`

```rust
assert_eq!(counts.counts.len(), 2);
```

The test verifies that *two distinct keys* exist but not:
- Which two keys they are (left AC / right GT as encoded)
- That each has weight 1.0

A future refactor that produces two identical keys (which would merge) or produces the
wrong keys would pass this test. The test should assert the specific key-value pairs.

---

### T2 — `encode_within_code` for `KmerSource::Read` is not directly unit-tested

**File:** `motifs_tests.rs`

The three direct tests of `encode_outside_code` and `encode_within_code` all use
`KmerSource::Reference`. The read-backed path (`spec.encode_kmer_bytes(&end.within_bases)`)
is exercised only indirectly through `count_fragment_in_window` tests. There is no test
that directly calls `encode_within_code` with `KmerSource::Read` and verifies the exact
encoded value.

---

### T3 — `encode_blacklist_validation_within_code` has no isolated test

**File:** `motifs.rs:449-476`

The blacklist-validation path — where read-source within encoding first validates the
reference to check for blacklisted bases before using the read sequence — has no unit
test. The only way to hit it is with a non-empty `blacklist_intervals` slice in the
context, which none of the current tests provide.

---

### T4 — `collect_end_motif_order` (sparse path) is not tested at all

**File:** `output.rs:31-39`, `output_tests.rs`

The four tests in `output_tests.rs` cover size guards and `build_all_end_motif_order`
(dense path only). `collect_end_motif_order` — the primary code path for the default
non-`--all-motifs` run — has zero direct tests. A sort-order regression or set-union
bug would go undetected.

---

### T5 — `stack_end_motif_counts` (dense matrix assembly) has no test

**File:** `write.rs:175-194`

The function that converts per-window sparse maps into a dense `ndarray` matrix is
never tested. A row/column transposition or motif-lookup regression would be invisible
until integration testing.

---

### T6 — `build_tile_payload` sort and serialization round-trip are not tested

**File:** `tiling.rs:88-116`, `tiling_tests.rs`

`tiling_tests.rs` has exactly one test (`merge_tile_payload_merges_counts_by_window_and_key`).
Missing:
- Merging into multiple distinct windows
- Empty payload
- `build_tile_payload` sort correctness (within-window and across-window)
- `serialize_tile_counts` / `deserialize_tile_counts` round-trip

---

### T7 — `decode_end_motif_counts` with `k_outside > 0` and `k_within = 0` not tested

**File:** `counting.rs` tests

All counting tests use `k_within > 0` with `outside_spec = None`. The case
`k_outside > 0, k_within = 0` (outside-only motif) is never exercised. Given the
asymmetric decode logic (`storage_order` differs for left/right ends), this is a
meaningful gap.

---

### T8 — `build_all_end_motif_order` with both `within_spec` and `outside_spec` not tested

**File:** `output_tests.rs`

Both output tests that exercise `build_all_end_motif_order` pass `None` for one of
the specs. No test verifies the combined case, where the label includes both an outside
and a within component, and collapse operates on the full `outside+within` concatenated
string.

---

### T9 — `merge_tile_payload` with multiple distinct windows not tested

**File:** `tiling_tests.rs`

The sole test merges two payloads for the *same* window index (7). It does not test
merging payloads that cover *different* windows (e.g., window 7 and window 12), which
exercises the `HashMap::entry` path for new keys, not just the additive path for
existing keys.

---

### T10 — `reference_motif_context_shares_code_table_when_within_and_outside_k_match` tests implementation detail

**File:** `motifs_tests.rs:235-246`

```rust
assert!(Arc::ptr_eq(within_codes, outside_codes));
```

This test asserts that the `Arc` pointers are identical — an internal memory
optimization. If the optimization is ever changed to a non-sharing strategy (still
semantically correct), the test fails. Tests should assert behavior (same codes at
same positions), not internal representation.

---

### T11 — `decode_end_motif_counts_collapses_reverse_complements_when_requested` test has potential misidentified pair

**File:** `counting.rs:355-383`

The test claims "GT" and "AC" are reverse complements. Verify:
`rev_complement("GT")` = reversed "TG", complement each → "AC". ✓
`rev_complement("AC")` = reversed "CA", complement each → "GT". ✓
Both canonicalize to "AC" (since lex("AC") < lex("GT")). ✓

The test is sound. Note however that `reverse_on_decode: false` is used for BOTH
keys — meaning both are treated as left ends. A complementary test for right-end keys
(where `reverse_on_decode: true` and the full decode includes a RC step) is missing.

---

### T12 — `decode_end_motif_counts_drops_motifs_with_n` constructs an N-containing key directly

**File:** `counting.rs:388-404`

```rust
within_code: within_spec.encode_kmer_bytes(b"AN"),
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

### D1 — `ClipStrategy::Raw` default in `EndsConfig::new()` is not documented

There is no doc comment or `//` note explaining why `new()` uses `Raw` while the CLI
uses `Aligned`. The discrepancy is invisible to programmatic users.

---

### D2 — `Endpoint` mode doc does not explain how it interacts with `CountOverlap` weight

**File:** `config_structs.rs:98-113`

The `--assign-by` documentation explains `count-overlap` as "count up the fraction of
fragment bases overlapping each window" but does not clarify that, in this mode, BOTH
ends of a fragment are counted in EVERY candidate window (not just the window the end
falls in). This is a significant semantic difference that users familiar with
`endpoint` mode may miss.

---

### D3 — `bins.bed` `overlap_perc` column meaning is undocumented in code or help text

No inline comment or help text explains what `overlap_perc` represents in the output
BED. Is it the fraction of the window covered by blacklisted regions? The fraction
of the window covered by the BAM? This is opaque.

---

### D4 — `--k-within = 0` and `--k-outside = 0` check is documented but the allowed combination (`k_within=0 OR k_outside=0`) is not

**File:** `ends.rs:74-76`

The guard allows one of the two to be zero. Using `k_within=0` means only outside
bases are counted; using `k_outside=0` means only within bases. The motif label format
becomes `"GT_"` or `"_AC"`. This is not documented in the help text for either
`--k-within` or `--k-outside`.

---

### D5 — `ClipStrategy::Drop` vs `--max-soft-clips 0` equivalence is not documented

Both cause soft-clipped ends to be skipped, but `Drop` uses `SkipEndDropAssignmentBoundary`
while `max_soft_clips = 0` also uses `SkipEndDropAssignmentBoundary`. For counting
purposes they are equivalent. This should be noted so users don't combine them
redundantly (or misunderstand that they're not subtly different).

---

### D6 — The `within_bases` orientation convention for `KmerSource::Read` on reverse reads is undocumented

**File:** `ends_fragment.rs:560-573`

`extract_right_within_bases` takes the last `k_within` bases of `read.seq`. This relies
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
intermediate steps could still misidentify what "outside" and "within" mean
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

### M7 — `EndsConfig::new()` silently accepts `k_within=0` and `k_outside=0` simultaneously

**File:** `config.rs:246-286`, `ends.rs:74-76`

The validation `k_within == 0 && k_outside == 0` happens only in `run()`, not in
`new()`. Constructing a config with both zero is possible programmatically without
an error until `run()` is called.

---

### M8 — `format_end_motif_label` called `format_end_motif_label` but just splits and rejoins with underscore

**File:** `counting.rs:185-196`

The function name is accurate but its docstring says "The full motif string is expected
to be oriented already and ordered as `outside || within`." This ordering assumption
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
format!("{outside}{within}")
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

For `k_within=1, k_outside=0, KmerSource::Read`:
- Left within = first 1 base of forward seq = `'A'` → label `_A`
- Right within = last 1 base of reverse seq = `'T'`, decoded as `RC('T') = 'A'` → label `_A`

Both endpoints produce the **same motif label**. Any integration test using this
fixture with `KmerSource::Read` cannot distinguish left from right endpoint by motif
label. A regression that swapped left/right endpoints would pass every such test as
long as the counts sum correctly.

The two CLI tests do not pass `--source-within`, so they use the Read default and both
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
  path in motif encoding, is only constructed in `drop_invalid_gc` (line 785) for GC
  purposes, not for motif encoding.

---

#### IT-I3 — `base_config` explicit `ClipStrategy::Aligned` override masks pre-fix B1

**File:** `test_ends_command.rs:61`

```rust
cfg.clip.clip_strategy = ClipStrategy::Aligned;
```

This was necessary to work around B1 (`EndsConfig::new()` defaulting to `Raw`). After
B1 is fixed, the override is redundant but harmless. More importantly: there is no
regression test asserting that `EndsConfig::new()` produces `ClipStrategy::Aligned`.
If B1 silently regresses (someone changes the default back to `Raw`), `base_config`
would still mask it and all tests would continue to pass.

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

#### IT-S1 — Indel filter tests check count, not identity of surviving ends

**File:** `test_ends_command.rs:875-878, 931-934`

`auto_indel_filter_keeps_indel_affected_read_backed_end_motifs`:
```rust
assert_eq!(matrix.shape(), &[1, 2]);  // 2 motifs
assert_eq!(motifs.len(), 2);
assert_eq!(matrix.sum(), 2.0);
```

`auto_indel_filter_skips_indel_affected_reference_backed_end_motifs`:
```rust
assert_eq!(matrix.shape(), &[1, 1]);  // 1 motif
assert_eq!(motifs.len(), 1);
assert_eq!(matrix.sum(), 1.0);
```

Neither test asserts **which end survived**. In the reference-backed case, the left
end (forward read, CIGAR `2M 1I 5M`, indel within the 4-base left footprint) should
be dropped, and the right end (reverse read, CIGAR `6M`, no indel) should survive.
But the test only checks `motifs.len() == 1`. A regression where the indel filter
accidentally kept the left end and dropped the right end would pass both assertions,
since either way there is one surviving motif. The correct assertion would compare
the surviving motif label against the expected right-end decode of the reference
sequence at positions 110–113.

---

#### IT-S2 — `all_requires_the_full_fragment_to_overlap_the_window` tests only rejection

**File:** `test_ends_command.rs:1202-1236`

Fragment `[10,20)` is tested against window `[10,19)`. Because `[10,19)` does not
fully contain `[10,20)`, the result is `sum = 0.0`. This confirms rejection. But there
is no test where the window fully contains the fragment (e.g., `[9,21)` or `[10,20)`)
to confirm that `All` mode actually counts when containment holds. A broken `All`
implementation that always rejects would pass all current tests.

---

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

Fragment B uses `KmerSource::Read` (no `--source-within` flag). Forward seq = AAAA,
reverse seq = TTTT. Left within = 'A' → `_A`. Right within = RC('T') = 'A' → `_A`.
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

- Left raw base: `seq[0]` of `b"TTACGTAA"` = `'T'`. Left within = `'T'` → label `_T`. ✓
- Right raw base: `seq[7]` = `'A'`. Right within = `'A'`, decoded as `RC('A') = 'T'` → label `_T`. ✓

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
- Outside-only label format: `"CGTAC_"` (outside bases then underscore for empty within). ✓

This test has the most rigorous mental derivation in the suite.

---

#### IT-C2 — `blacklist_masking_still_skips_read_backed_within_motifs` derivation note

**File:** `test_ends_command.rs:430-432`

```
// - left read-backed motif = "_A"
// - right read-backed motif = reverse-complement("A") = "_T"
```

Read seq = `b"ACGA"`, so:
- Left within = `seq[0] = 'A'` → `_A`. ✓
- Right within = `seq[3] = 'A'`, decoded as `RC(b'A') = b'T'` → `_T`. ✓

The comment writes "reverse-complement("A") = '_T'", which applies RC to the *base*
`'A'` to get `'T'`, then forms the label `_T`. This is correct but could be read as
applying RC to the motif string `"A"` (which also gives `"T"` for a 1-mer). The
description is unambiguous for k=1 but would be wrong for k>1 if interpreted as
motif-level RC. Minor documentation risk for future test authors.

---

### 8.4 Missing Coverage

#### IT-M1 — `WindowMotifAssigner::All` acceptance case is not tested

**File:** No test file

The only `All`-mode integration test (`all_requires_the_full_fragment_to_overlap_the_window`,
line 1203) tests rejection: window `[10,19)` is 1 bp smaller than fragment `[10,20)`,
so nothing is counted. There is no test where the window is large enough to fully
contain the fragment (e.g., `[10,20)` or `[9,21)`), so the acceptance branch of the
`All` predicate has no integration-level coverage. A bug where `All` always rejects
would pass all current tests.

---

#### IT-M2 — No integration test for `k_within > 0 AND k_outside > 0` simultaneously

**File:** No test file

Every integration test sets either `k_within=1, k_outside=0` (the vast majority) or
`k_within=0, k_outside=5` (the fallback test). The combined case — where both within
and outside bases are present and the label is `"outside_within"` — is never exercised
end-to-end. The label format, the encoding ordering (outside concatenated before within
for the left end), and the collapse-complement behavior on full concatenated labels are
all untested in integration.

---

#### IT-M3 — `collapse_complement` with odd `k_within + k_outside ≥ 3` is not tested

**File:** No test file

The `collapse_complement_merges_reverse_complement_equivalent_end_motifs` test (line
1390) uses `k_within=1, k_outside=0`. Total k=1 (odd). For k=1, the middle-base
heuristic and lexicographic canonicalization coincide: RC('G') = 'C' and 'C' < 'G',
so both methods return 'C'. The two approaches can diverge for odd k≥3 (issue L15).
No integration test exercises `collapse_complement=true` with k_within+k_outside ≥ 3
and odd, so the L15 heuristic discrepancy is entirely untested.

---

#### IT-M4 — No multi-chromosome integration test

**File:** All tests use `base_chromosomes(&["chr1"])`

All 56 tests restrict processing to a single chromosome. The chromosome-iteration
logic in `run()`, `chr_offsets_for_threads`, and the per-chromosome scaling lookup
(`scaling_chr` in `ends.rs:510`) are never exercised with more than one chromosome.
A regression in chromosome switching or off-by-one in tile assignment across chromosome
boundaries would be invisible.

---

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

#### IT-M6 — No test for the B6 error path (incomplete scaling coverage)

**File:** No test file

`scaling_factors_weight_each_counted_end_motif` (line 470) provides a scaling file
covering `[0, 256)` — the entire chromosome. There is no test where the scaling file
has a gap and a fragment falls in an uncovered region. The B6 issue (opaque error
message) is therefore not regression-tested.

---

#### IT-M7 — No regression test for B1 (ClipStrategy default after fix)

**File:** No test file

B1 is now fixed, but there is no test asserting `EndsConfig::new().clip.clip_strategy
== ClipStrategy::Aligned`. Since `base_config` in the test suite explicitly sets
`Aligned`, a future regression to `Raw` would not be caught by the integration tests.

---

#### IT-M8 — No test for proportion value precision in settings JSON (L14)

**File:** No test file

`settings_json_records_clip_and_window_assignment_semantics` (line 2143) does not use
`WindowMotifAssigner::Proportion`. The L14 risk — that `proportion=0.2` could
serialize as `"proportion=0.19999999999999998"` — has no test.

---

#### IT-M9 — `global_mode` tested with only one fragment

**File:** `test_ends_command.rs:1489`

`global_mode_counts_both_end_motifs_in_one_output_row` uses a single fragment.
"Global mode" (no windowing) is the default when no `--by-size` or `--by-bed` is
specified. With a single fragment the output must be a 1×N matrix. There is no test
confirming that two fragments produce correctly accumulated counts in the same single
output row.

---

#### IT-M10 — No test for the `bins.bed` 4th-column value (C7)

**File:** No test file

`windowed_runs_write_bins_bed_with_the_selected_windows` (line 2228) asserts:
```rust
assert!(rows[0].starts_with("chr1\t10\t11\t"));
```
The 4th column (overlap_perc from C7) is present in the output but never tested for
correctness. Its value is undocumented and unvalidated. A regression changing the 4th
column from a float to an integer, or from overlap percentage to some other statistic,
would pass all current tests.

---

#### IT-M11 — No integration test for `k_within` or `k_outside` > 1 with KmerSource::Reference

**File:** No test file

All within-only tests use `k_within=1`. The fallback test uses `k_outside=5` with
`KmerSource::Read` (not `KmerSource::Reference`). No test exercises reference-backed
within-bases with `k_within > 1`. The `get_reference_code` logic for multi-base
within reference lookup — including sentinel_n handling for any N in the k-mer — has
no integration coverage.

---

### 8.5 Structural Issues

#### IT-ST1 — Settings JSON tests use `contains()` only — cannot detect malformed JSON

**File:** `test_ends_command.rs:2143-2180, 2268-2293`

```rust
assert!(settings.contains("\"k_within\": 1"));
assert!(settings.contains("\"clip_strategy\": \"aligned\""));
```

These assertions confirm field presence but not JSON validity. A regression that
introduced a trailing comma (breaking JSON syntax), duplicated a field, or omitted a
field not currently asserted would go undetected. The tests should at minimum parse the
JSON with `serde_json::from_str` to confirm it is valid, and ideally assert on the
complete field set.

---

#### IT-ST2 — `by_size_windowing_writes_bins_bed` only checks "non-empty" and "chr1\t" prefix

**File:** `test_ends_command.rs:2506-2508`

```rust
assert!(!bins_bed.trim().is_empty());
assert!(bins_bed.lines().all(|row| row.starts_with("chr1\t")));
```

This does not verify the number of windows generated, their sizes, or whether the
window size matches the requested `by_size = Some(20)`. The test passes for any
non-empty chr1 BED output, including incorrect window boundaries.

---

#### IT-ST3 — Output prefix existence tests don't verify absence of un-prefixed files

**File:** `test_ends_command.rs:2404-2439, 2542-2571`

The prefix tests verify that `sampleA.end_motifs.sparse.npz` etc. exist. They do not
verify that `ends.end_motifs.sparse.npz` (or any un-prefixed variant) was NOT created.
A regression that wrote both the prefixed and un-prefixed copies would pass.

---

#### IT-ST4 — `both_kmer_sizes_zero_is_rejected` does not validate error message

**File:** `test_ends_command.rs:2125-2141`

```rust
assert!(run(&cfg).is_err());
```

Only checks that an error is returned, not what message it carries. This is inconsistent
with `unpaired_mode_rejects_require_proper_pair` (line 2314), which asserts the error
contains `"--require-proper-pair cannot be used with --reads-are-fragments"`. A future
change that returned the wrong error (e.g., a downstream panic rather than a validation
error) would not be caught.

---

#### IT-ST5 — CLI tests don't exercise `--source-within reference`

**File:** `test_ends_command.rs:263-340, 345-424`

Both CLI tests omit `--source-within`. The CLI defaults to `KmerSource::Read` (no
reference required). No CLI test exercises the full pipeline through the binary with
`--source-within reference`, meaning the CLI argument parsing for that flag is not
integration-tested at the binary level. Only the programmatic API tests use
`KmerSource::Reference`.

---

#### IT-ST6 — All fixture fragments use the same start/length geometry

**File:** Most tests use `simple_paired_fragment_bam("...", 10, 10, 4)`

With rare exceptions (edge tests, cross-tile), every standard fragment is at position
10 with aligned length 10. Bugs that only manifest for fragments at position 0,
fragments whose k_outside window crosses the chromosome start, or fragments with
lengths that interact with tile-size modulo arithmetic would not be caught by the
standard geometry.

---

### 8.6 Test Suite Positive Observations

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
| [√] B1 | Bug | High | config.rs, config_structs.rs |
| [√] B2 | Bug | Medium | counters.rs, ends.rs |
| M11 | Style | Low | ends.rs |
| B4 | Bug | Low | ends.rs |
| [√] D8 | Docs | Low | ends.rs |
| [√] B6 | Bug | Medium | ends.rs |
| [√] L1 | Logic | Low | ends.rs |
| L2 | Logic | Medium | ends.rs |
| [√] L3 | Logic | Low | ends.rs |
| [√] L4 | Logic | Low | ends.rs, motifs.rs |
| L5 | Logic | Low | motifs.rs |
| L6 | Logic | Low | motifs.rs |
| M12 | Style | Low | counting.rs |
| [√] L8 | Logic | Low | write.rs |
| L9 | Logic | Low | output.rs |
| L10 | Logic | Medium | config_structs.rs |
| L11 | Logic | Low | motifs_tests.rs |
| L12 | Logic | Low | ends_fragment.rs |
| L13 | Logic | Low | ends_fragment.rs |
| [√] L14 | Logic | Medium | write.rs |
| L15 | Logic | High | base.rs, counting.rs |
| C1 | Config | Low | config.rs |
| C3 | Config | Medium | config_structs.rs |
| C4 | Config | Low | write.rs |
| C5 | Config | Medium | write.rs |
| C6 | Config | Low | write.rs |
| C7 | Config | Medium | ends.rs |
| C8 | Config | Medium | config_structs.rs |
| [√] T1 | Tests | Medium | motifs_tests.rs |
| [√] T2 | Tests | Medium | motifs_tests.rs |
| [√] T3 | Tests | High | motifs.rs |
| [√] T4 | Tests | High | output_tests.rs |
| [√] T5 | Tests | High | write.rs |
| T6 | Tests | Medium | tiling_tests.rs |
| T7 | Tests | Medium | counting.rs |
| T8 | Tests | Medium | output_tests.rs |
| T9 | Tests | Low | tiling_tests.rs |
| T10 | Tests | Low | motifs_tests.rs |
| T11 | Tests | Low | counting.rs |
| T12 | Tests | Low | counting.rs |
| P1 | Perf | Low | motifs.rs |
| P2 | Perf | Low | ends.rs |
| P3 | Perf | Low | ends.rs |
| P4 | Perf | Low | tiling.rs |
| D1 | Docs | High | config.rs |
| D2 | Docs | Medium | config_structs.rs |
| D3 | Docs | Low | ends.rs |
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
| IT-I3 | Tests/Infra | Low | test_ends_command.rs |
| IT-I4 | Tests/Infra | Low | fixtures/mod.rs |
| IT-S1 | Tests | Medium | test_ends_command.rs:829,883 |
| IT-S2 | Tests | Medium | test_ends_command.rs:1202 |
| IT-S3 | Tests | Low | test_ends_command.rs:1084 |
| IT-S4 | Tests | Low | test_ends_command.rs:1153 |
| IT-S6 | Tests | Low | test_ends_command.rs:422 |
| IT-M1 | Tests | High | (no test) |
| IT-M2 | Tests | High | (no test) |
| IT-M3 | Tests | Medium | (no test) |
| IT-M4 | Tests | Medium | (no test) |
| IT-M5 | Tests | Medium | (no test) |
| IT-M6 | Tests | Low | (no test) |
| IT-M7 | Tests | Low | (no test) |
| IT-M8 | Tests | Low | (no test) |
| IT-M9 | Tests | Low | (no test) |
| IT-M10 | Tests | Low | (no test) |
| IT-M11 | Tests | Medium | (no test) |
| IT-ST1 | Tests | Low | test_ends_command.rs |
| IT-ST2 | Tests | Low | test_ends_command.rs:2506 |
| IT-ST3 | Tests | Low | test_ends_command.rs:2426 |
| IT-ST4 | Tests | Low | test_ends_command.rs:2125 |
| IT-ST5 | Tests | Medium | test_ends_command.rs:263,345 |
| IT-ST6 | Tests | Low | test_ends_command.rs |
