# midpoints Spec

`cfdna midpoints` builds grouped midpoint profiles around fixed-width genomic sites. The public output is dense, while tile internals are sparse.

## Input Contract

- `--intervals` is required and must be BED-like with columns `chromosome`, `start`, `end`, and `group_name`.
- The public CLI expects sorted intervals. The current loader sorts selected windows internally, but group index assignment still follows first observed group name in the input file.
- All selected intervals must have the same length.
- Sites with the same group name collapse into one group profile.
- Group indices are assigned by first observed group name during grouped BED loading.
- Chromosome filtering may remove rows, but the command must fail if no selected grouped windows remain.

## Fragment And Midpoint Geometry

- Paired fragments span `[forward.pos, reverse.reference_end)`.
- Unpaired read-as-fragment mode spans `[read.pos, read.reference_end)`.
- Length bins are half-open and also define the accepted fragment length range.
- Even-length fragment midpoints use deterministic coordinate-derived random rounding. Duplicate fragments with the same coordinates choose the same midpoint base.
- Tile ownership is by midpoint position in the tile core, not by fragment start.
- Blacklist filtering is fragment-level and happens before midpoint ownership is checked in the current implementation.

## Counting

- Each accepted fragment contributes to every grouped site window containing its midpoint.
- Position is zero-based relative to the site's interval start.
- The dense public array has shape `(group, length_bin, position)`.
- The flat storage order is group-major, then length-bin, then position.
- Fragment length lookup must use the shared `LengthAxis`, not a manual bin search.
- Weight per contribution is GC weight times optional scaling weight.
- Scaling is averaged over the full aligned fragment and applied to every selected midpoint window.

## GC, Scaling, And Blacklists

- GC correction can come from `--gc-file` plus `--ref-2bit`, or from a two-byte BAM AUX `--gc-tag`.
- Invalid GC weights skip fragments by default. `--neutralize-invalid-gc` keeps them with weight 1.0.
- File-based GC correction reads reference sequence for the full tile fetch span.
- Scaling TSVs must fully cover every selected chromosome and pass GC-mode compatibility checks.
- Blacklist strategy uses the aligned fragment span and the shared blacklist strategies: `any`, `all`, `midpoint`, or `proportion=<threshold>`.
- `midpoint` blacklist strategy checks the single central base for odd fragments and either central base for even fragments. It does not use the randomized counted-midpoint tie break.

## Tiling And Merge

- The command first selects grouped site windows overlapping a tile core, then narrows BAM fetch to the extreme selected site span plus maximum-fragment-length halo, clamped to tile fetch.
- Tile workers count into `SparseProfileGroupsCounts`.
- Non-empty tile partials are written as sparse `.npz` archives with arrays:
  - `idx.npy`: sorted flat dense indices as `u64`
  - `data.npy`: `f32` counts
  - `shape.npy`: `[group, length_bin, position]`
- Empty sparse tiles do not write partial files.
- Merge allocates one final dense output buffer and adds sparse partials in parallel into chunk locks.
- Sparse partials must validate shape, sorted indices, equal `idx`/`data` lengths, platform index fit, and destination bounds before merging.

## Output Contract

- Main output is `<prefix>.midpoint_profiles.npy` with shape `(group, length_bin, position)`.
- Group index output is `<prefix>.group_index.tsv`.
- With the `plotters` feature, selected group indices can emit quick QC line plots and length-bin heatmaps.
- Run statistics include counted fragments, blacklist exclusions, and GC failure summaries when relevant.

## Open Notes

! Warning: Projection from paired sites to midpoint-centered or strand-aware profiles remains future work. Keep projection plans outside finalized specs until implemented.

! Warning: File-based GC prefix construction still happens before no-window pruning for sparse targeted runs. This is a performance issue, not a correctness mismatch.
