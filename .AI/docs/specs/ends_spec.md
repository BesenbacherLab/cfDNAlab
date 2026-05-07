# ends Spec

`cfdna ends` counts fragment end motifs. It is a motif-counting command, but its core risk is geometry: clipped boundaries, outside reference lookup, motif-level blacklist validation, and window assignment must stay aligned.

## Motif Definition

- Each end motif is labeled as `<outside>_<inside>`.
- `k_inside` and `k_outside` may be zero independently, but both cannot be zero for useful output.
- The left end is decoded as `outside || inside`.
- The right end is decoded from stored `inside || outside`, then reverse-complemented so final labels are oriented from the fragment end inward in reference 5'->3' orientation.
- Motifs containing `N` or sentinel invalid bases are dropped during final decoding.
- `--source-inside read` uses sequenced inside bases.
- `--source-inside reference` uses reference bases and requires `--ref-2bit`.
- Outside bases always come from the reference and require `--ref-2bit` when `k_outside > 0`.
- `--collapse-complement` is experimental and hidden unless built with `ends_experimental`. It canonicalizes same-orientation complements after the full motif is joined.

## Clipping

- `skip` drops soft-clipped motifs at the affected end. This is the default.
- `aligned` ignores soft-clipped bases and uses aligned boundaries.
- `raw-aligned-boundary` uses raw read bases for inside sequence but keeps aligned genomic boundaries for outside lookup, window assignment, and motif-level blacklist validation.
- `raw-shifted-boundary` uses raw read bases and shifts the fragment-end boundary outward by the soft-clipped length.
- Hard-clipped fragments are always discarded.
- Raw clipping modes are supported only with `--source-inside read`.
- File-based GC correction and genomic scaling still use aligned reference span even when motif assignment uses shifted raw boundaries.
- For raw-shifted `count-overlap`, clipped-only window contributions use the nearest aligned reference base for scaling.

## Indels And Base Quality

- `--indel-filter auto` resolves to allow indel-affected motifs for read-backed inside bases and to skip affected ends for reference-backed inside bases.
- `skip-affected-end` drops only the end whose motif overlaps an indel.
- `skip-affected-fragment` drops the full fragment when either end is affected.
- `--bq-filter` requires `k_inside > 0` and `--source-inside read`.
- Base-quality filter grammar is `<agg> in <scope> <op> <threshold>`.
- Aggregation is `min`, `mean`, or `max`.
- Scope is `end` or `fragment`.
- Operators are `>=`, `>`, `<=`, and `<`.
- Repeated filters are conjunctive. End filters drop individual ends. Fragment filters drop the full fragment.

## Blacklists

- Fragment-level blacklist filtering uses the assignment geometry implied by the clip strategy.
- Motif-level blacklist validation always drops motifs whose reference-addressable motif bases overlap masked reference sequence.
- Read-backed inside motifs still use read bases, but when a blacklist is active the reference-addressable inside portion is validated against masked reference.
- In `raw-aligned-boundary`, clipped-only inside bases have no reference coordinate and are not checked by motif-level blacklist validation.
- `--ref-2bit` is required when a blacklist is active and motif validation needs reference bases.

## Window Assignment

- Global mode counts all kept motifs into one row.
- Fixed-size, BED, and grouped BED modes use `WindowMotifAssigner`.
- `endpoint` assigns each end independently by its own endpoint position. The right endpoint is `boundary_pos - 1` in half-open coordinates.
- `count-overlap` counts fragment overlap fraction in each selected window.
- `any`, `all`, `midpoint`, and `proportion=<threshold>` select windows by fragment assignment interval and count full end motif mass in selected rows.
- Even-length midpoint assignment uses deterministic coordinate-derived random rounding.
- Grouped BED windows map output rows to group indices and collapse counts by group.

## Weighting

- Each counted end receives: assignment overlap mass, times GC weight, times scaling weight.
- GC weight is fragment-level and reused for all counted ends and windows from that fragment.
- Scaling for non-`count-overlap` modes is averaged over the full aligned fragment and reused for selected windows.
- Scaling for `count-overlap` is averaged over the aligned overlap span for each output row.
- Scaling TSVs must fully cover selected chromosomes and pass GC-mode compatibility checks.

## Tiling And Temporary Counts

- Tile ownership is by fragment aligned start in the tile core. This prevents double counting across fetch halos.
- BED fetch narrowing uses candidate-window extent with a halo at least as large as maximum fragment length.
- Motif reference preload currently uses the full tile fetch band plus outside/raw padding, independent of narrowed BAM fetch.
- Tile workers write sorted sparse bincode payloads keyed by global window or group id.
- Reduction merges returned tile payload paths directly. It does not scan temporary directories for work.
- Per-fragment statistics count `counted_fragments` once when at least one end motif contributed, and `counted_motifs` as the number of distinct ends that contributed.

## Output Contract

- Sparse output writes `<prefix>.end_motifs.sparse.npz` plus `<prefix>.end_motifs.txt`.
- Dense output with `--all-motifs` writes `<prefix>.end_motifs.npy` plus `<prefix>.end_motifs.txt`.
- Dense output enumerates every possible motif label and is guarded by `CFDNALAB_ENDS_MAX_DENSE_OUTPUT_BYTES`, default 5 GiB.
- Motif labels are sorted deterministically.
- Settings sidecar is `<prefix>.end_motif_settings.json` and records motif lengths, inside source, clip strategy, window assignment, indel filter, effective indel filter, base-quality filters, and experimental complement collapse when enabled.

## Open Notes

! Warning: The settings sidecar intentionally does not yet record every normalization input, such as GC/scaling file paths and dense-vs-sparse format. Do not treat it as a complete provenance record.

! Warning: Motif reference preload is still full-tile plus motif padding even when sparse BED fetch narrowing skips most aligned reads. This is a performance issue, not a correctness mismatch.
