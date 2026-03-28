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

- one encoded `within` part
- one encoded `outside` part
- one orientation flag that tells output whether the combined motif must be reverse-complemented on reconstruction

These are counted in a dedicated struct and only turned into full motif strings when writing output.

This keeps counting fast and compact while preserving the exact motif identity.

### 4. Outside and within should use the same radix-5 encoding family

For consistency and compactness:

- `outside` reference context should use the same positional k-mer encoding style as `fragment-kmers`
- `within` sequence should use the same radix-5 encoding, but encoded directly from the resolved end sequence rather than from a chromosome-wide positional table

### 5. Complement collapsing happens on the combined full motif, not per half

If `collapse_complement` is enabled, it must be applied to the full joined motif after `within` and `outside` are combined and oriented correctly.

Do not canonicalize the two halves independently.

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
    pub within_bases: Vec<u8>,
}
```

Notes:

- `boundary_pos` is the assignment boundary used by the kept end
- for `aligned`, `boundary_pos` is the aligned boundary after terminal clipping is excluded
- for `raw`, `boundary_pos` may move outward beyond the aligned span while the counted `within_bases` still come from the raw terminal sequence
- if an end is skipped, no `ResolvedFragmentEnd` exists for it

This shape is intentionally minimal. If later implementation shows that another small field is needed for counters or validation, it can be added.

## Per-read collection design

### `EndReadInfo`

As in the other fragment collectors, per-read collection should first extract a compact read summary and then combine one or two summaries into `FragmentWithEnds`.

The per-read summary needs enough information to:

- orient the read relative to the fragment
- inspect terminal clipping
- build the resolved within-fragment sequence for the end
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
- resolve `within_bases` from the first aligned bases inside the fragment
- clipped terminal bases are excluded from the counted within-fragment sequence

This is the default clip strategy.

### `raw`

Collector behavior:

- keep the end
- resolve `within_bases` from the observed read bases at the fragment end, including soft-clipped bases
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

These are derived from the first `k_within` aligned bases used by the selected clip strategy. The collector does not need to store indel offsets or other intermediate state.

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
  - if the end motif has indels in its aligned `k_within` footprint, return
    `SkipEndKeepAssignmentBoundary`

- `IndelMotifFilterPolicy::Auto`
  - in `Reference` source mode, treat indel-affected motifs the same as `SkipAffectedEnd`
  - in `Read` source mode, keep the end

- `IndelMotifFilterPolicy::SkipAffectedFragment`
  - if either end motif has indels in its aligned `k_within` footprint, return `DropFragment`

## Required tests

These tests are required before calling the collector stable enough.

### Edge clipping parsing

- left soft clip only, e.g. `10S90M`
- right soft clip only, e.g. `90M10S`
- hard clip only at left, e.g. `5H95M`
- hard clip only at right, e.g. `95M5H`
- combined hard+soft edge patterns, e.g. `5H10S90M` and `90M10S5H`

### Clip strategy behavior

- `aligned` excludes terminal soft-clipped bases from `within_bases`
- `raw` includes terminal soft-clipped bases in `within_bases`
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
- indels outside the aligned `k_within` footprint do not trigger skipping

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

The within-fragment side should use the same radix-5 encoding family, but encoded directly from the resolved end sequence.

This means:

- do not build chromosome-wide within-position tables
- encode only the short `within_bases` slice already resolved by the collector

This keeps the representation uniform without pretending that within-fragment sequence is a reference-position track.

## Count key design

Each counted end motif should be stored internally in a dedicated key struct:

```rust
pub struct EncodedEndMotifKey {
    pub within_code: u64,
    pub outside_code: u64,
    pub reverse_complement_on_decode: bool,
}
```

The window-level counter should use this struct as the key.

The important point is that counting is done on encoded motif parts plus decode-time orientation, not on strings.

## Decoding and output

At reduction/output time:

1. decode `within_code` using the `k_within` spec
2. decode `outside_code` using the `k_outside` spec
3. join them into the full motif string in storage order
4. if `reverse_complement_on_decode` is true, reverse-complement the joined full motif
5. if `collapse_complement` is enabled, canonicalize the oriented full motif
6. aggregate duplicate full motifs if canonicalization merges them

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
5. Add direct encoding for resolved `within` sequences
6. Add positional reference encoding for `outside`
7. Count encoded `(within_code, outside_code)` tuples
8. Decode and combine motifs only at final output

## Explicit non-goals for this step

- implementing `fill-with-ref`
- prematurely reusing `fragment-kmers` payloads
- carrying raw clipping ambiguity far downstream
- decoding full motif strings in the hot counting loop
- canonicalizing motif halves independently
- reverse-complementing every fragment end in the hot loop just to normalize orientation before counting
