# Ends Motifs File Spec

This is a future design note for adding `--motifs-file` behavior to `cfdna ends`.
The goal is to let users restrict counting to a known motif subset, and optionally
count directly into user-defined motif groups, so large combined end motifs can be
handled without expanding memory usage to the full motif universe.

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

The output target order in file mode should follow file order. For grouped files,
group order should be first-seen order.

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

The file path should precompute a lookup from encoded counted state to output
target id:

```text
EncodedEndMotifKey -> target_id
```

The key must include all fields:

```rust
EncodedEndMotifKey {
    inside_code,
    outside_code,
    reverse_on_decode,
}
```

Including `reverse_on_decode` is required. Reverse-complement-related motifs can
belong to different user groups, so they must not be accidentally collapsed or
matched through decoded string equivalence.

For each parsed motif label, precompute the exact encoded left-end key and the
exact encoded right-end key that decode to that final oriented label. Both keys
map to the same file-defined target for that row. Do not generate additional
reverse-complement aliases unless those exact labels are present as separate
rows in the file.

## Counting Shape

Keep the no-file branch simple and fast. A small mode enum is enough:

```rust
enum EndMotifSelection {
    None,
    Selected(SelectedEndMotifLookup),
}
```

At the point where a validated `EncodedEndMotifKey` would normally be counted,
branch once:

```rust
match selection {
    EndMotifSelection::None => {
        encoded_window_counts.incr_weighted(encoded_key, weight);
    }
    EndMotifSelection::Selected(lookup) => {
        if let Some(target_id) = lookup.target_for(encoded_key) {
            selected_window_counts.incr_weighted(target_id, weight);
        }
    }
}
```

The selected path should count into compact numeric target ids instead of storing
the full encoded motif keys. This is what gives grouped motif files their memory
benefit.

Avoid making the existing `EndMotifCounts` generic unless that becomes clearly
cleaner during implementation. A parallel selected counter is likely easier to
read:

```rust
struct SelectedEndMotifCounts {
    counts: FxHashMap<u32, f64>,
}
```

The helper that stores one observation should return whether a motif was actually
counted, so `CountedEndFlags` remains correct when a motifs file filters out an
otherwise valid end motif.

## Validation Rules

Validate the motifs file before opening BAM-derived outputs.

Fail on:

- Invalid number of columns.
- Empty motif labels or empty group names.
- Invalid bases.
- Motif lengths that do not match `--k-outside` and `--k-inside`.
- Invalid underscore use for the current inside/outside lengths.
- Motifs containing `N`.
- Duplicate motif rows.
- The same exact encoded key mapping to multiple target ids.
- `--motifs-file` combined with `--collapse-complement`.

If the same exact encoded key appears twice with the same target id, prefer
failing anyway. Duplicate rows usually indicate a bad input file.
