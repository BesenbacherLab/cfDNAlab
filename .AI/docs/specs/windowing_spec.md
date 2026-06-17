# windowing Spec

Windowing selects output rows. Fetch narrowing selects aligned reads. These are related but separate decisions, and future command work must keep them separate.

## Window Modes

- `WindowSpec::Global`: one output row.
- `WindowSpec::Size(bp)`: contiguous fixed windows per chromosome from 0 to chromosome length. Final window is clipped.
- `WindowSpec::Bed(path)`: ordinary BED rows with original row indices.
- `DistributionWindowSpec::GroupedBed(path)`: grouped BED rows where row identity becomes `group_idx`.
- Distribution grouped BED behaves like ordinary BED for fetch/window geometry, but downstream output identity differs.

## BED Loading

- BED coordinates are 0-based half-open `[start, end)`.
- Empty lines, comments, `track`, and `browser` directives are skipped.
- Valid windows increment original BED index even if filtered out by chromosome selection or a filter function.
- Ordinary BED stores `(start, end, original_idx)`.
- Grouped BED stores `(start, end, group_idx)` and requires a fourth column.
- Group indices are assigned by first observed group name.
- Per-chromosome window lists are start-sorted after loading.
- Commands must fail early when selected BED/grouped BED windows are empty after chromosome filtering.

## Window IDs

- Global mode maps every chromosome to row 0.
- Fixed-size mode maps chromosome-local bin indices through per-chromosome offsets.
- BED mode uses embedded original BED indices and offset 0.
- Grouped BED mode uses embedded group indices.
- `build_bin_info` sorts BED metadata by original index before writing.

## Assignment Semantics

Shared `WindowAssigner`:

- `count-overlap`: count the fraction of fragment bases overlapping each window.
- `any`: select windows with any fragment-position overlap.
- `all`: select windows where the fragment positions are effectively fully inside the window.
- `midpoint`: select windows containing the deterministic fragment midpoint.
- `proportion=<threshold>`: select windows meeting a fragment-position overlap fraction.

`ends` uses `WindowMotifAssigner`, which adds `endpoint`:

- `endpoint`: assign each motif by its own fragment-end coordinate.
- Other modes are fragment-centric and then count kept end motifs into selected rows.

For proportion modes, the denominator is fragment positions, not window positions. A small window can fully cover itself and still cover only a small proportion of a long fragment.

## Fetch Policy

Use `fetch_span_for_tile` rather than ad hoc BED-to-fetch logic.

- Global mode returns the tile fetch span, clamped to chromosome length.
- Fixed-size mode derives the fixed-window extent touching the tile core, then clamps with halo.
- BED mode must choose one of `CoreOverlap`, `CandidateWindowExtent`, or `KeepTileFetch`.
- `window_derived_fetch_extent_for_core_overlap` must use true tile-core overlaps, even when a cached candidate span is wider.
- `window_derived_fetch_extent_for_candidates` must use the candidate span directly and must not reapply core-overlap filtering.

## Metadata Files

- Ordinary window metadata uses `chrom`, `start`, `end`, and `blacklisted_fraction`.
- Grouped metadata always includes `group_idx` and `group_name`.
- Grouped blacklist fractions are interval-width weighted across all intervals in the group.
- Overlapping grouped intervals contribute separately to both numerator and denominator.
- Group names in TSV output must be sanitized for tabs and newlines.

## Implementation Invariants

- Window slices passed to streaming helpers must be start-sorted.
- Window offsets and original/group indices are public output contracts. Do not reindex BED rows in reducers unless the command explicitly writes a new metadata file with that mapping.
- Keep candidate-window selection, fetch narrowing, fragment ownership, and row identity as separate concepts in code and docs.
