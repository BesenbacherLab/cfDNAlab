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
