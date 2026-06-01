# Ends Motifs File Spec

This is a tracking note for the partially implemented `--motifs-file` behavior
in `cfdna ends`. The goal is to let users restrict counting to a known motif
subset, and optionally count directly into user-defined motif groups, so large
combined end motifs can be handled without expanding memory usage to the full
motif universe.

## Current Gap

The selected counting/output path is wired, but selected-subspace encoding is
not wired into the command yet.

Currently implemented:

- The parser accepts one-column motif files and two-column grouped motif files.
- Counting can route selected observations into compact numeric target columns.
- Tile payloads, reduction, post-processing, and Zarr writing have selected
  paths.
- `--all-motifs` with `--motifs-file` uses the file-defined target axis rather
  than the full motif universe.

Not yet implemented:

- The motifs-file path still builds ordinary `KmerSpec` values through
  `build_optional_kmer_spec`.
- That keeps the ordinary per-side k limit and prevents `--motifs-file` from
  counting k-mers larger than the radix-5 full-space representation supports.
- `SubspaceKmerSpec` exists in `kmer_codec.rs`, but it is not used by
  `cfdna ends`.

Until this is fixed, the CLI help promise about much larger motifs is not fully
true.

## Implementation Plan

1. Parse the motifs file into validated motif rows before building encoded
   lookup keys. Keep target assignment and group first-seen order as currently
   implemented.

2. Build selected side subspaces from the parsed rows:
   - one optional `SubspaceKmerSpec` for inside bases when `k_inside > 0`
   - one optional `SubspaceKmerSpec` for outside bases when `k_outside > 0`
   - no subspace for a disabled side; the disabled side should still encode as
     the existing zero code.

3. Replace the selected lookup key with a key whose side codes come from the
   selected subspaces, not from ordinary radix-5 `KmerSpec` codes. Keep
   `reverse_on_decode` in the key.

4. Encode parser-derived left and right observable states through the selected
   subspace specs:
   - left key uses the motif as written
   - right key uses reverse-complemented inside/outside halves and
     `reverse_on_decode = true`
   - both keys map to the target assigned by the motif row.

5. Update selected counting so it encodes observed motif halves through the
   selected subspace specs and skips immediately when either enabled side returns
   the missing sentinel. Only construct the combined selected key after both
   enabled sides are selected.

6. Update reference precomputation for selected motifs. Tile context should
   build one selected-code array per enabled side from `SubspaceKmerSpec`.
   Radix-backed subspaces can reuse ordinary precomputation and remap codes;
   byte-backed subspaces scan the reference slice and look up normalized bytes.

7. Keep the no-file path unchanged. It should continue using ordinary
   `KmerSpec`, ordinary `EncodedEndMotifKey`, existing decoding, and the current
   k limit.

8. Add regression tests:
   - `--motifs-file` with `k_inside > 27` succeeds on a tiny selected axis.
   - unselected large-k observations are skipped before target insertion.
   - right-end reverse-complement state maps to the correct target/group.
   - no-file large k still fails on the ordinary path.
   - `--motifs-file + --collapse-complement` is rejected.

9. After the implementation is real, distill current behavior into
   `.AI/docs/specs/ends_spec.md` and remove or rewrite this future note so it no
   longer describes stale design-only structures.

## Motif Identity

The motifs file defines combined end motifs, not separate outside and inside
motif sets.

Rows use the same public motif identity as the `ends` output:

```text
<outside>_<inside>
```

Each side must match `--k-outside` and `--k-inside`. When only one side has
non-zero length, the underscore may be omitted by the user, but the parsed motif
should be normalized internally to the command's public combined motif identity.

The file labels are interpreted after the normal `ends` orientation rule: each
motif is in fragment-end-inward 5'->3' orientation. The file must not be treated
as a list of reverse-complement-equivalent motifs.

## File Columns

The file is tab-separated.

```text
motif
motif	group
```

With one column, each motif is its own output target.

With two columns, motifs are counted directly into the named group. This grouping
is the intended replacement for automatic complement collapsing when users want
custom motif collapses.

Group names should be restricted to ASCII letters, numbers, `.`, `_`, and `-`.
Do not allow whitespace in group names.

The output target order in file mode should follow file order. For grouped files,
group order should be first-seen order.

## Selected Targets

The parsed motifs file should produce one output target axis, regardless of
whether the file is grouped.

```rust
struct SelectedEndMotifLookup {
    labels: Vec<String>,
    column_kind: EndMotifColumnKind,
    inside_spec: Option<SubspaceKmerSpec>,
    outside_spec: Option<SubspaceKmerSpec>,
    lookup: FxHashMap<SelectedEndMotifKey, TargetId>,
}
```

For an ungrouped file, each row creates one motif target:

```text
AT_CG -> target 0, label "AT_CG"
GC_TA -> target 1, label "GC_TA"
```

For a grouped file, rows create or reuse group targets:

```text
AT_CG	cpg_like -> target 0, label "cpg_like"
GC_TA	cpg_like -> target 0, label "cpg_like"
TT_AA	other    -> target 1, label "other"
```

The combined selected lookup maps each counted end state to the final output
target id. It should not preserve motif-level identities in grouped mode:

```text
left encoded key for AT_CG  -> target 0
right encoded key for AT_CG -> target 0
left encoded key for GC_TA  -> target 0
right encoded key for GC_TA -> target 0
left encoded key for TT_AA  -> target 1
right encoded key for TT_AA -> target 1
```

The target labels are the output axis labels. In grouped mode, motif labels from
the file are validation and lookup inputs, not output columns.

`column_kind` records whether those labels are motifs or motif groups. It should
be `Motif` for one-column files and `MotifGroup` for two-column files.

The target id is an index into `SelectedEndMotifLookup::labels`. Its storage
type should be chosen from the target-axis size rather than hardcoded to `u32`.
Small target spaces should not pay for wider indices than needed.

## Subspace K-mer Encoding

The motifs-file path should not use the full-space `build_kmer_specs` limit as
its core representation. The file defines a selected k-mer subspace separately
for the inside side and the outside side:

```rust
struct SubspaceKmerSpec {
    k: usize,
    code_dtype: CodeDType,
    encoding: SubspaceKmerEncoding,
    missing_code: CodeValue,
}

enum SubspaceKmerEncoding {
    Radix5 {
        full_spec: KmerSpec,
        selected_by_radix_code: FxHashMap<u64, CodeValue>,
    },
    Bytes {
        selected_by_bytes: FxHashMap<Box<[u8]>, CodeValue>,
    },
}
```

For k values that fit the ordinary radix-5 representation, the subspace spec can
reuse `KmerSpec` to encode observed kmers quickly and then remap full-space
codes to compact subspace codes. For k values that do not fit radix-5 `u64`, the
subspace spec should use byte-slice lookup.

The compact subspace code is side-specific. A combined selected end motif key is
then:

```rust
struct SelectedEndMotifKey {
    inside_code: CodeValue,
    outside_code: CodeValue,
    reverse_on_decode: bool,
}
```

The code value type should be adaptive. If the inside or outside subspace has
few kmers, `u8` or `u16` is enough. Larger subspaces may need `u32` or `u64`.
This mirrors the existing `KmerCodes` width selection. Do not hardcode selected
side codes or target ids to `u32` if a narrower type is sufficient or a wider
type is required.

Each subspace spec has a missing sentinel. Invalid, N-containing, out-of-bounds,
or unselected kmers map to that sentinel, allowing selected counting to skip
before constructing a combined target lookup key.

Reference precomputation should produce one selected-code array per side, not
separate left and right arrays. A k-mer code at a reference position can be used
by either fragment direction depending on the requested genomic start and the
`reverse_on_decode` state carried by the final combined key. This matches the
ordinary reference-code precomputation model.

## Interaction With Existing Options

`--motifs-file` and `--all-motifs` are compatible.

Without `--motifs-file`, `--all-motifs` keeps its existing meaning: dense output
over the full possible motif universe.

With `--motifs-file`, `--all-motifs` means dense output over the file-defined
target axis. For ungrouped files this is the listed motifs. For grouped files
this is the listed groups.

`--motifs-file` should not be allowed with `--collapse-complement`. The file can
already express arbitrary grouping, and automatic complement collapsing would
make validation and interpretation unnecessarily ambiguous.

## Lookup Semantics

Filtering must happen before inserting into per-window count maps. Filtering only
at decode or output time would not prevent memory growth during counting.

The normal no-file path should stay close to the current implementation: extract
and count `EncodedEndMotifKey` values directly, then decode later for output.

The file path should precompute side subspaces and a combined lookup from
selected side codes to output target id:

```text
SelectedEndMotifKey -> target_id
```

The key must include all fields:

```rust
SelectedEndMotifKey {
    inside_code,
    outside_code,
    reverse_on_decode,
}
```

Including `reverse_on_decode` is required. Reverse-complement-related motifs can
belong to different user groups, so they must not be accidentally collapsed or
matched through decoded string equivalence.

For each parsed motif label, precompute the exact selected-code left-end key and
the exact selected-code right-end key that decode to that final oriented label.
Both keys map to the same file-defined target for that row. Do not generate
additional reverse-complement aliases unless those exact labels are present as
separate rows in the file.

## Counting Shape

Keep the no-file branch simple and fast. The command only needs an optional
selected-motif lookup:

```rust
Option<SelectedEndMotifLookup>
```

At the point where a validated `EncodedEndMotifKey` would normally be counted,
the selected path should instead look up side codes through the subspace specs
and skip when either side returns the missing sentinel:

```rust
match selected_motifs {
    None => {
        encoded_window_counts.incr_weighted(encoded_key, weight);
    }
    Some(lookup) => {
        let Some(inside_code) = lookup.inside_code_for_end(...) else { return; };
        let Some(outside_code) = lookup.outside_code_for_end(...) else { return; };
        let key = SelectedEndMotifKey { inside_code, outside_code, reverse_on_decode };
        if let Some(target_id) = lookup.target_for(key) {
            selected_window_counts.incr_weighted(target_id, weight);
        }
    }
}
```

The selected path should count into compact numeric target ids instead of storing
the full encoded motif keys. This is what gives grouped motif files their memory
benefit.

For reference-backed sides, build tile-local selected-code arrays from the
subspace specs:

```text
reference position -> selected side code or missing sentinel
```

For radix-backed subspaces, this can reuse ordinary radix-5 precomputation and
remap full-space codes into selected side codes. For byte-backed subspaces, scan
the reference slice and look up each k-mer in the byte map. This path may be
slower for very large k, but avoids full-universe memory growth.

For read-backed inside motifs, use the inside subspace spec directly on
`ResolvedFragmentEnd::inside_bases`. Blacklist validation remains a separate
reference-span check and should not be responsible for producing the selected
inside code.

Avoid making the existing `EndMotifCounts` generic unless that becomes clearly
cleaner during implementation. A selected count map is enough:

```rust
FxHashMap<TargetId, f64>
```

The helper that stores one observation should return whether a motif was actually
counted, so `CountedEndFlags` remains correct when a motifs file filters out an
otherwise valid end motif.

## Output Format

New end-motif outputs should use schema version 2 and declare the column-axis
kind in the root metadata:

```json
"motif_axis_kind": "motif"
```

or:

```json
"motif_axis_kind": "motif_group"
```

Existing schema version 1 stores are ordinary motif-axis stores without
`motif_axis_kind`; downstream readers may keep accepting those as old motif
outputs.

Ungrouped output targets are ordinary motifs. This includes no-file output and
one-column motifs files. They should keep the existing `motif_ascii` metadata
layout, but still be written as schema version 2 so older downstream packages
fail clearly and users update.

The count arrays should keep their existing shapes and dimension names:

```text
counts[row, motif]
sparse/{row,motif,count,shape,sparse_dimension}
```

Grouped output targets are motif groups. For motif-group output, `motif_index`
is still the numeric count-column axis, but its labels are stored directly in
JSON metadata:

```json
{
  "label_field": "motif_group",
  "labels": ["group.one", "group-two"]
}
```

Do not write `motif_byte` or `motif_ascii` for motif-group output. Group labels
can be variable width and should not be zero-padded into the motif ASCII matrix.

Downstream readers should keep the schema distinction but not duplicate the
public selector API. `motifs()` and motif selectors refer to the count-column
axis for both axis kinds. When `motif_axis_kind` is `motif_group`, those labels
are user-defined group names rather than concrete DNA motifs.

## Validation Rules

Validate the motifs file before opening BAM-derived outputs.

Fail on:

- Invalid number of columns.
- Empty motif labels or empty group names.
- Group names containing anything other than ASCII letters, numbers, `.`, `_`, and `-`.
- Invalid bases.
- Motif lengths that do not match `--k-outside` and `--k-inside`.
- Invalid underscore use for the current inside/outside lengths.
- Motifs containing `N`.
- Duplicate motif rows.
- The same exact encoded key mapping to multiple target ids.
- `--motifs-file` combined with `--collapse-complement`.

If the same exact encoded key appears twice with the same target id, prefer
failing anyway. Duplicate rows usually indicate a bad input file.
