# Ends implementation spec

Date: 2026-03-28

## Scope

This spec covers the first implementation shape for `cfdna ends`.

It is intentionally narrow. It only fixes the internal representation and counting model needed to get the command working in a clean way. It does not try to settle every future feature now.

The main goals are:

- define the fragment payload shape that `ends` should stream from the fragment iterator
- decide where clipping is resolved
- define how motifs are represented internally for counting
- defer `fill-with-ref` until the simpler path is working

## Decision summary

### 1. `ends` gets its own fragment payload type

`ends` should follow the existing codebase pattern:

- command-specific fragment payload in `src/shared/fragment/`
- matching per-read info type
- matching `collect_fragment_with_...` logic
- matching iterator entry point in `src/shared/fragment_iterator.rs`

The command should not try to reuse `FragmentWithIndelCounts` or `FragmentWithKmerSegments` as-is.

Working name:

- `FragmentWithEnds`

### 2. Clipping is resolved during fragment collection

The collector should apply the selected clip strategy and only pass downstream the information that remains valid after that decision.

This keeps downstream logic simple and fits the repo pattern where fragment payloads already bake in command-relevant preprocessing.

This means:

- `aligned`, `raw`, and `drop` are handled in fragment collection
- `fill-with-ref` is not implemented yet and should explicitly use `unimplemented!()`

### 3. Internal counting uses encoded motif halves plus orientation

Each counted end motif should be represented internally as:

- one encoded `inside` part
- one encoded `outside` part
- one orientation flag that tells output whether the combined motif must be reverse-complemented on reconstruction

These are counted in a dedicated struct and only turned into full motif strings when writing output.

This keeps counting fast and compact while preserving the exact motif identity.

### 4. Outside and inside should use the same radix-5 encoding family

For consistency and compactness:

- `outside` reference context should use the same positional k-mer encoding style as `fragment-kmers`
- `inside` sequence should use the same radix-5 encoding, but encoded directly from the resolved end sequence rather than from a chromosome-wide positional table

Implementation consequence:

- `shared/kmers/kmer_codec.rs` should remain the one codec implementation
- `ends` should not introduce its own parallel encoding logic
- the shared codec needs one small additional helper for direct encoding of a single byte slice, because the current public API is centered on chromosome-wide positional precomputation

Working addition:

- `KmerSpec::encode_kmer_bytes(&self, seq: &[u8]) -> u64`

This should:

- encode one exact k-mer-sized byte slice using the same radix-5 representation already used elsewhere
- return the existing sentinel values when the slice contains `N` or is otherwise not a valid full k-mer
- let `ends` encode `inside_bases` directly without duplicating codec logic

### 4b. `k_inside` or `k_outside` may be zero, but not both

`ends` should allow one side of the motif to be empty:

- `k_inside = 0`, `k_outside > 0`
- `k_inside > 0`, `k_outside = 0`

But it should reject:

- `k_inside = 0` and `k_outside = 0`

because that would define an empty motif and make the output meaningless.

Implementation rule:

- only build `KmerSpec` values and positional code tables for sides with `k > 0`
- when one side has `k = 0`, always store code `0` for that side in the internal key
- on decode/output, treat that side as the empty string rather than as a real k-mer

This keeps the count key compact and avoids `Option<u64>` in the hot path.

### 5. Complement collapsing happens on the combined full motif, not per half

If `collapse_complement` is enabled, it must be applied to the full joined motif after `inside` and `outside` are combined and oriented correctly.

Do not canonicalize the two halves independently.

For `ends`, the public contract is that the decoded motif string is already in final
biological `outside || inside` order before any collapsing happens. That means:

- `decode_full_motif` is the only step that may reverse the joined sequence
- `collapse_complement` must compare the decoded motif against its direct
  same-orientation complement, not against its reverse complement
- the canonical representative is the lexicographically smaller of
  `{motif, complement(motif)}`

Why this matters:

- `revcomp(outside || inside) = revcomp(inside) || revcomp(outside)`
- so reverse-complement canonicalization would swap the two conceptual halves
- the final split by `k_outside` would then put the wrong bases into the
  `outside` and `inside` slots

This is not a per-half rule. `outside || inside` is one combined end motif identity.
The full decoded motif is compared as one string against its same-orientation
complement, and only then split into the public `<outside>_<inside>` label.

Worked example:

- `k_outside = 1`, `k_inside = 2`
- left end decodes to `GTA`
- matching right-end storage decodes to `CAT`
- `complement("GTA") = "CAT"`
- `CAT < GTA`, so the canonical full motif is `CAT`
- final label is `C_AT`

At no point should the canonicalization step return `revcomp("GTA")`, because that
would be `TAC` and would no longer respect the `outside || inside` contract.

### 6. `ends` needs its own encoded count key

`shared/kmers/kmer_codec::Kmer` is not the right internal count key for `ends`.

It only models one k-mer plus orientation, while `ends` needs:

- one `inside` code
- one `outside` code
- one decode-time orientation flag

So `ends` should define its own dedicated key type in command-local code.

Working shape:

```rust
pub struct EncodedEndMotifKey {
    pub inside_code: u64,
    pub outside_code: u64,
    pub reverse_on_decode: bool,
}
```

Notes:

- this is an internal counting key, not a user-facing output type
- it should be `Eq + Hash + Clone + Copy`
- it should stay minimal and not accumulate unrelated fragment metadata
- complement collapsing still happens later on the full reconstructed motif, not on this key directly
- when `k_inside == 0`, `inside_code` must always be `0`
- when `k_outside == 0`, `outside_code` must always be `0`

## Fragment payload design

### `FragmentWithEnds`

The fragment payload should carry:

- fragment-level information needed for:
  - length filtering
  - blacklist filtering
  - window assignment
  - GC correction
  - genomic scaling
- resolved end-level information needed for motif extraction and counting

The payload should support both paired and unpaired input.

The command-level semantics are:

- paired input contributes two fragment ends from the aligned fragment span `forward.pos -> reverse.reference_end`
- unpaired input contributes two fragment ends from the aligned read span `pos -> reference_end`

So both modes use aligned fragment boundaries by default. They differ only in whether the fragment span comes from a read pair or from one read treated as a full fragment.

So the payload should represent both fragment ends for both input modes without forcing downstream code to branch on raw read layouts.

Working shape:

```rust
pub struct FragmentWithEnds {
    pub tid: i32,
    pub interval: Interval<u32>,
    pub gc_tag: GcTagValue,
    pub left_end: Option<ResolvedFragmentEnd>,
    pub right_end: Option<ResolvedFragmentEnd>,
}
```

This shape is preferred over a small vector because:

- left and right end identity stays explicit
- paired and unpaired input use the same downstream representation
- later orientation handling does not need an extra per-end role tag

### `ResolvedFragmentEnd`

Each resolved end should carry only what downstream counting needs after clip handling has already been applied.

Working shape:

```rust
pub struct ResolvedFragmentEnd {
    pub boundary_pos: u32,
    pub inside_bases: Vec<u8>,
}
```

Notes:

- `boundary_pos` is the assignment boundary used by the kept end
- for `aligned`, `boundary_pos` is the aligned boundary after terminal clipping is excluded
- for `raw`, `boundary_pos` may move outward beyond the aligned span while the counted `inside_bases` still come from the raw terminal sequence
- if an end is skipped, no `ResolvedFragmentEnd` exists for it

This shape is intentionally minimal. If later implementation shows that another small field is needed for counters or validation, it can be added.

## Counting model

### Sparse per-window counting

`ends` should follow `fragment-kmers`, not `lengths`, for the count container shape.

`lengths` is still a good model for:

- tiling
- GC correction
- blacklist and scaling plumbing
- reduction over tile outputs

But its dense one-dimensional count arrays are the wrong shape for motifs.

So `ends` should use sparse per-window maps keyed by `EncodedEndMotifKey`.

## Output conventions

### Dense vs sparse final output

Final output should only be dense when `--all-motifs` is enabled.

Reason:

- dense output is useful for small motif spaces where users want the same motif columns across samples
- when `--all-motifs` is disabled, the natural output is sparse because only observed motifs need to be written

Implementation rule:

- `all_motifs = true`:
  - enumerate the full motif universe
  - write a dense `.npy` matrix
  - estimate the dense output size up front against a configurable output-size budget and fail early if it would exceed that budget
- `all_motifs = false`:
  - collect only observed motifs
  - write a sparse COO `.npz` matrix plus the matching motif label file

The dense-output budget should be based on actual output size, not only on `k`.

Working threshold:

- default to roughly 4-5 GiB, with an environment-variable override for users who want to allow larger dense outputs

### Motif labels

The final motif labels should be written as:

- `<outside>_<inside>`

Examples:

- `GG_ATC`
- `_AC` when `k_outside = 0`
- `GG_` when `k_inside = 0`

Implementation rule:

- first decode and orient the full motif in the existing `outside || inside` order
- then, if `collapse_complement` is enabled, compare that full decoded motif against
  its same-orientation complement and keep the lexicographically smaller of the two
- then split by `k_outside`
- then format directly as the user-facing `<outside>_<inside>` label

This keeps the current motif orientation and collapse semantics intact while making the displayed order match the cut-centered motif orientation.

### Settings sidecar

The settings sidecar should be renamed from the inherited `lengths` filename.

Working filename:

- `end_motif_settings.json`

It should include the settings needed to interpret the output, especially:

- `k_inside`
- `k_outside`
- `source_inside`
- `clip_strategy`
- `max_soft_clips`
- `indel_filter`
- `window_assignment`
- `collapse_complement`
- `reads_are_fragments`
- fragment-length filter settings

### Fragment-length filtering

Fragment-length filters remain defined on the aligned fragment interval, not on the assignment interval.

Working shape:

```rust
type EndMotifCounts = FxHashMap<EncodedEndMotifKey, f64>;
type WindowEndMotifCounts = FxHashMap<u64, EndMotifCounts>;
```

Meaning:

- one tile accumulates sparse counts per global window id
- each window stores only the motifs actually seen in that tile
- reduction merges sparse maps across tiles before final decoding/output

This is the same overall strategy already used by `fragment-kmers`.

### Decode only at the end

The hot path should count encoded keys only.

During reduction/output:

- decode `inside_code`
- decode `outside_code`
- join them into the full motif
- apply `reverse_on_decode`
- then apply optional same-orientation complement collapse on the already oriented
  full motif

This keeps the streaming path small and fast.

### `ends` should get a small dedicated counting module

Add a small command-local counting module, for example:

- `src/commands/ends/counting.rs`

This module should own:

- `EncodedEndMotifKey`
- sparse per-window count types
- small helpers for incrementing, merging, and later decoding counts

It should not duplicate generic radix-5 logic from `shared/kmers/`.

## Per-read collection design

### `EndReadInfo`

As in the other fragment collectors, per-read collection should first extract a compact read summary and then combine one or two summaries into `FragmentWithEnds`.

The per-read summary needs enough information to:

- orient the read relative to the fragment
- inspect terminal clipping
- build the resolved inside-fragment sequence for the end
- support unpaired and paired input consistently

It should include:

- `tid`
- aligned interval
- strand
- terminal clip lengths relevant to the fragment end
- oriented read-end sequence needed for `raw` and `aligned`
- GC tag if used

The exact per-read shape should be chosen to match the existing collector style in:

- `segment_fragment.rs`
- `segment_kmer_fragment.rs`
- `indel_counting_fragment.rs`

### SAM strand convention that this implementation must follow

This needs to be stated explicitly because it is easy to get wrong.

For mapped reads, the BAM/SAM record sequence already follows the reference orientation stored in the record:

- mapped segments are represented on the forward genomic strand
- if a read is aligned to the reverse strand, `SEQ` is already reverse-complemented relative to the original unmapped read
- for reverse-strand alignments, strand-sensitive fields are also reversed consistently with that stored sequence

For `ends`, BAM storage orientation is not yet the final motif orientation we want to count.

The BAM record stores mapped reverse-strand reads reverse-complemented so they follow the reference genome. But `ends` wants end-specific motif orientation, meaning the sequence should run from the fragment end inward.

So for `ends`:

- use `is_reverse` to decide which side of the stored `SEQ` corresponds to the fragment end of interest
- then reverse-complement when needed to convert the extracted end sequence from BAM/reference-oriented storage into end-specific motif orientation

In other words:

- BAM `SEQ` on reverse-strand alignments is already complemented to follow the reference genome
- `ends` must often reverse-complement the extracted end sequence again so it is oriented from the fragment end inward

This is not "fixing BAM". It is converting from BAM storage orientation to the motif orientation that this command counts.

## Clip strategy boundary

### Default fragment-end interpretation

By default, `ends` should trust the aligner and use aligned fragment ends.

This applies in both input modes:

- paired fragments use the aligned span `forward.pos -> reverse.reference_end`
- unpaired `--reads-are-fragments` input uses the aligned span `pos -> reference_end`

So clipped terminal bases are not part of the fragment by default. Users opt into using clipped terminal sequence only by selecting `raw`.

### `aligned`

Collector behavior:

- keep the end
- use the aligned fragment boundary as the counted end
- resolve `inside_bases` from the first aligned bases inside the fragment
- clipped terminal bases are excluded from the counted inside-fragment sequence

This is the default clip strategy.

### `raw`

Collector behavior:

- keep the end
- resolve `inside_bases` from the observed read bases at the fragment end, including soft-clipped bases
- if terminal soft clipping is present, move the counted fragment boundary outward beyond the aligned span by the clipped length
- clipping is allowed

### `drop`

Collector behavior:

- if the relevant fragment end is clipped, omit that end
- paired fragments may still keep the opposite end
- unpaired fragments may therefore yield zero ends and be skipped entirely

### `fill-with-ref`

Current decision:

- do not implement yet
- use `unimplemented!()` in the collector path

Reason:

- the simpler path should work first
- this mode has the weakest semantics and should be re-evaluated later

## Indels and blacklist interaction

This spec does not redefine those semantics. It assumes the current API direction already chosen for `ends`:

- blacklist has two layers:
  - fragment-level exclusion by `blacklist_strategy`
  - motif-level invalidation always on
- indel handling for motifs is controlled by `IndelMotifFilterPolicy`

The collector should therefore resolve the end-level filtering outcomes directly.

For indels, the collector only needs motif-level booleans:

- `left_motif_has_indels`
- `right_motif_has_indels`

These are derived from the first `k_inside` aligned bases used by the selected clip strategy. The collector does not need to store indel offsets or other intermediate state.

### Blacklist handling must be two-stage

`ends` must follow the same two-stage blacklist idea as `fragment-kmers`, but adapted to end motifs.

Stage 1:

- apply fragment-level exclusion using `is_blacklisted(...)` on the aligned fragment interval
- this uses the configured `blacklist_strategy`
- if the fragment is excluded here, stop before any end motif work

Stage 2:

- every end motif must also be checked against a blacklist-masked reference representation
- this stage is always on, independent of `blacklist_strategy`
- if the genomic span used by an end motif touches a blacklisted base, that end motif is skipped

The fast path must not use per-end overlap queries. Instead, it should reuse the `fragment-kmers` idea:

- load the tile reference sequence when blacklist-aware motif validation is needed
- mask blacklisted bases in that tile slice before precomputing k-mer codes
- use the existing radix-5 sentinel behavior to detect invalid spans

Implementation rule:

- a blacklist-masked motif span should produce the same invalid signal as any other ambiguous base span
- in practice this means the masked-reference code at the lookup point should decode to the existing `N` sentinel and the end should be skipped

This rule applies to both motif halves:

- `outside` always validates against the masked reference span for that outside lookup
- `inside` also validates against the masked reference span for the inside lookup, even when the actual inside sequence is taken from the read

So for `source_inside = read`:

- the actual `inside_code` still comes from `inside_bases`
- but the corresponding genomic inside span must first be validated against the masked reference
- if that masked reference span produces the `N` sentinel, skip the end instead of counting the read-derived motif

This keeps blacklist semantics genomic without forcing expensive overlap checks on every end.

### End-specific lookup coordinates for blacklist validation

Blacklist validation must use the same genomic lookup starts as the reference-backed motif encoding.

For a kept end with assignment boundary `boundary_pos`:

- left `outside` span starts at `boundary_pos - k_outside`
- left `inside` span starts at `boundary_pos`
- right `outside` span starts at `boundary_pos`
- right `inside` span starts at `boundary_pos - k_inside`

The right-end offsets are easy to get wrong and must be tested explicitly.

### Exact fallback path must preserve blacklist semantics

If a requested reference-backed motif span falls outside the currently preloaded tile slice and the implementation falls back to an exact reference fetch, that fallback must preserve blacklist masking semantics.

In other words:

- exact fallback must not read raw reference sequence and bypass the blacklist mask
- the fallback result must behave exactly like the masked tile-precomputed path for the same genomic span

## Windowing and run-shape

### Reuse the existing window helpers

For window indexing and metadata, `ends` should follow the current `fragment-kmers` helpers rather than copying another custom windowing path.

The useful existing pieces are:

- `WindowContext`
- `compute_window_offsets(...)`
- `build_bin_info(...)`

These already handle:

- `--by-size`
- `--by-bed`
- `--global`
- chromosome-local to global window index mapping
- output row ordering

So the first implementation should reuse them directly, even if they later move to a more neutral shared location.

### Overall runner shape

For the top-level `run()` implementation:

- follow `lengths` for tile processing, GC correction, blacklist/scaling plumbing, and reduction flow
- follow `fragment-kmers` for sparse counting and final decode/write steps

In other words:

- `lengths` is the right model for tile orchestration
- `fragment-kmers` is the right model for motif counting representation

This split is intentional and should be preserved while implementing `ends.rs`.

### Tile ownership and downstream windows

The sparse tile payload design must preserve the same fragment ownership invariant as the existing tiled commands.

Rule:

- a fragment is counted only in the tile whose core contains the fragment start

In practice:

- every tile fetches a wider halo region so it can see fragments near its boundaries
- a fragment visible in multiple fetch halos must still be counted only once
- the ownership check is the fragment-start-in-core rule

This ownership rule is separate from window assignment.

Once a tile owns a fragment, that fragment may still contribute to windows outside the tile core if the chosen window-assignment rule allows it.

So the implementation must also preserve:

- window lookup spans that extend beyond the core far enough to include all windows that an owned fragment can legitimately hit

This means:

- sparse tile payloads do not need the old dense "cross" file format
- but they still must include counts for downstream windows beyond the tile core
- later reduction should merge sparse counts by global `original_idx`
- duplicate counting should be prevented by tile ownership, not by the reducer

This is important and should not be "simplified away" later. A fragment can be owned by one tile while still contributing to windows that begin in a later tile.

## Immediate next steps

The implementation should proceed in this order:

1. add the direct byte-slice radix-5 encoder to `shared/kmers/kmer_codec.rs`
2. add the `ends` counting module with `EncodedEndMotifKey` and sparse per-window counts
3. wire `ends.rs` to use the shared window helpers and sparse counting path
4. only then add the tile-level motif extraction loop and output writing

## Assignment boundary rules

The command now has two distinct geometries:

- aligned fragment interval geometry
  - used for fragment length and GC correction
  - always stays on the original aligned interval

- assignment boundary geometry
  - used for end-based placement when an end is actually kept for assignment
  - may differ from the aligned boundary under `raw`

This distinction is important and must stay explicit in the collector.

### End-resolution outcomes

One resolved end should return one of these semantic outcomes:

- `DropFragment`
  - abort the whole fragment
  - used for policies such as `skip-affected-fragment`

- `SkipEndKeepAssignmentBoundary { assignment_boundary_pos }`
  - skip the motif for this end
  - still let this end contribute a clip-adjusted or otherwise meaningful assignment boundary

- `SkipEndDropAssignmentBoundary`
  - skip the motif for this end
  - do not let this end redefine assignment geometry
  - callers should fall back to the aligned boundary

- `KeepEnd { assignment_boundary_pos, end }`
  - keep both the motif and the assignment boundary

### Current rules for choosing the outcome

- `aligned`
  - if the end can be built, return `KeepEnd` with aligned assignment boundary
  - if too few sequence bases remain, return `SkipEndDropAssignmentBoundary`

- `raw`
  - if the end can be built, return `KeepEnd` with outward-shifted assignment boundary
  - if too few sequence bases remain after using raw sequence, return `SkipEndDropAssignmentBoundary`

- `drop`
  - if the end is soft-clipped, return `SkipEndDropAssignmentBoundary`
  - otherwise behave like `aligned`

- `max_soft_clips`
  - when exceeded, return `SkipEndDropAssignmentBoundary`

- `IndelMotifFilterPolicy::SkipAffectedEnd`
  - if the end motif has indels in its aligned `k_inside` footprint, return
    `SkipEndKeepAssignmentBoundary`

- `IndelMotifFilterPolicy::Auto`
  - in `Reference` source mode, treat indel-affected motifs the same as `SkipAffectedEnd`
  - in `Read` source mode, keep the end

- `IndelMotifFilterPolicy::SkipAffectedFragment`
  - if either end motif has indels in its aligned `k_inside` footprint, return `DropFragment`

## Required tests

These tests are required before calling the collector stable enough.

### Edge clipping parsing

- left soft clip only, e.g. `10S90M`
- right soft clip only, e.g. `90M10S`
- hard clip only at left, e.g. `5H95M`
- hard clip only at right, e.g. `95M5H`
- combined hard+soft edge patterns, e.g. `5H10S90M` and `90M10S5H`

### Clip strategy behavior

- `aligned` excludes terminal soft-clipped bases from `inside_bases`
- `raw` includes terminal soft-clipped bases in `inside_bases`
- `raw` moves the assignment boundary outward by the soft-clipped length
- `drop` skips a soft-clipped end
- hard-clipped reads are always discarded

### Boundary semantics

- fragment `interval` stays on aligned boundaries for length and GC purposes
- kept `raw` ends use shifted assignment boundaries
- skipped ends from `drop` or `max_soft_clips` fall back to aligned assignment boundaries
- skipped ends from indel filtering keep assignment boundaries when policy says to keep them

### Indel filtering

- `Auto` + `Read` source keeps indel-affected ends
- `Auto` + `Reference` source skips indel-affected ends but keeps assignment boundaries
- `SkipAffectedEnd` skips only the affected end
- `SkipAffectedFragment` drops the whole fragment
- indels outside the aligned `k_inside` footprint do not trigger skipping

### Paired vs unpaired consistency

- paired fragments expose left and right ends from `forward.pos -> reverse.reference_end`
- unpaired fragments expose left and right ends from `pos -> reference_end`
- equivalent clipped layouts in paired and unpaired input produce the same end-level outcomes

## Reference context encoding

### Outside bases

The outside-fragment side should use reference-sequence positional encoding like `fragment-kmers`.

Use:

- `KmerSpec`
- `build_left_aligned_codes_per_k(...)`
- tile-relative code lookup from the fetched reference slice

This should provide encoded `k_outside` motifs at genomic positions without repeated string work.

### Within bases

The inside-fragment side should use the same radix-5 encoding family, but encoded directly from the resolved end sequence.

This means:

- do not build chromosome-wide inside-position tables
- encode only the short `inside_bases` slice already resolved by the collector

This keeps the representation uniform without pretending that inside-fragment sequence is a reference-position track.

## Count key design

Each counted end motif should be stored internally in a dedicated key struct:

```rust
pub struct EncodedEndMotifKey {
    pub inside_code: u64,
    pub outside_code: u64,
    pub reverse_complement_on_decode: bool,
}
```

The window-level counter should use this struct as the key.

The important point is that counting is done on encoded motif parts plus decode-time orientation, not on strings.

## Decoding and output

At reduction/output time:

1. decode `inside_code` using the `k_inside` spec
2. decode `outside_code` using the `k_outside` spec
3. join them into the full motif string in storage order
4. if `reverse_complement_on_decode` is true, reverse-complement the joined full motif
   so the result is now in final biological `outside || inside` order
5. if `collapse_complement` is enabled, canonicalize that already oriented full motif by
   taking the lexicographically smaller of `{motif, complement(motif)}`
6. aggregate duplicate full motifs if canonicalization merges them
7. split the canonical full motif at `k_outside` and write the final
   `<outside>_<inside>` label

This produces the final motif columns written to disk.

## Window assignment model

This spec assumes the current `ends` API direction:

- `endpoint` assigns each motif by its own fragment-end position
- the fragment-based assignment modes assign windows at the fragment level and then count the fragment's end motifs in those windows

This means the fragment payload must preserve:

- fragment interval for fragment-based assignment
- per-end boundary position for endpoint assignment

## Suggested implementation order

1. Add `FragmentWithEnds` and the matching per-read summary in `src/shared/fragment/`
2. Add `collect_fragment_with_ends(...)`
3. Add iterator support in `src/shared/fragment_iterator.rs`
4. Make the collector support:
   - `raw`
   - `drop`
   - `align-start`
   - `fill-with-ref => unimplemented!()`
5. Add direct encoding for resolved `inside` sequences
6. Add positional reference encoding for `outside`
7. Count encoded `(inside_code, outside_code)` tuples
8. Decode and combine motifs only at final output

## Explicit non-goals for this step

- implementing `fill-with-ref`
- prematurely reusing `fragment-kmers` payloads
- carrying raw clipping ambiguity far downstream
- decoding full motif strings in the hot counting loop
- canonicalizing motif halves independently
- reverse-complementing every fragment end in the hot loop just to normalize orientation before counting
