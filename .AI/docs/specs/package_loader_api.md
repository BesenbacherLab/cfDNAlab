# Package Loader API Spec

This spec describes the current public loader APIs for the Python and R helper
packages. The helper packages read files already produced by the Rust `cfdna`
CLI. They do not install, wrap, or reimplement the CLI.

## Shared Rules

- Expose selectors only for axes that exist on the loaded object. Global outputs
  do not take window or group selectors.
- Vector-capable selectors use plural names: `groups`, `group_idxs`,
  `window_idxs`, `motifs`, `motif_idxs`, `with_lengths`, and
  `length_bin_idxs`.
- Length-range selectors use `with_length_range` and take one half-open
  `[start, end)` bp range.
- Scalar inputs are accepted for plural selectors.
- Name selectors and index selectors for the same axis are mutually exclusive.
- Duplicate requested names or indices are errors.
- Name lookup checks uniqueness for the requested name. A duplicate unrelated
  name in the metadata must not poison a lookup for a unique name.
- Python indices follow ordinary Python conventions. R-facing indices and
  public `*_idx` columns are one-based.
- Public data frames must not expose internal zero-based R fields such as
  `*_idx0` or `*_index0`.
- Window metadata uses `chrom`, `start`, and `end` for public genomic ranges.

## Length Counts

Python:

```python
lengths = cfl.read_lengths("sample.length_counts.tsv.zst")

lengths.length_bins()
lengths.length_bin_idx(167)
lengths.counts_array(
    with_lengths=None,
    with_length_range=None,
    length_bin_idxs=None,
) # shape is always (output_row, length_bin)
lengths.window_metadata()    # windowed only
lengths.group_metadata()     # grouped only

lengths.data_frame(
    with_lengths=None,
    with_length_range=None,
    length_bin_idxs=None,
    value="count",
    denominator="all_bins",
    keep_wide=False,
)  # global
lengths.data_frame(
    window_idxs=None,
    with_lengths=None,
    with_length_range=None,
    length_bin_idxs=None,
    value="count",
    denominator="all_bins",
    keep_wide=False,
    max_blacklisted_fraction=1.0,
)  # windowed
lengths.data_frame(
    groups=None,
    group_idxs=None,
    with_lengths=None,
    with_length_range=None,
    length_bin_idxs=None,
    value="count",
    denominator="all_bins",
    keep_wide=False,
    max_blacklisted_fraction=1.0,
)  # grouped
```

R:

```r
lengths <- read_lengths("sample.length_counts.tsv.zst")

length_bins(lengths)
length_bin_idx(lengths, 167)
length_counts_vector(lengths) # global only
length_counts_matrix(lengths) # all modes

length_data_frame(
  lengths,
  with_lengths = NULL,
  with_length_range = NULL,
  length_bin_idxs = NULL,
  value = "count",
  denominator = "all_bins",
  keep_wide = FALSE
) # global
length_data_frame(
  lengths,
  window_idxs = NULL,
  with_lengths = NULL,
  with_length_range = NULL,
  length_bin_idxs = NULL,
  value = "count",
  denominator = "all_bins",
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0
) # windowed
length_data_frame(
  lengths,
  groups = NULL,
  group_idxs = NULL,
  with_lengths = NULL,
  with_length_range = NULL,
  length_bin_idxs = NULL,
  value = "count",
  denominator = "all_bins",
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0
) # grouped
```

Length data frames and arrays use the selected output rows and length bins.
Long output has one row per selected output unit and length bin. Wide output
has one row per selected output unit and one value column per length bin.

Length-bin selectors mean:

- `with_lengths`: bins containing exact fragment lengths. Multiple requested
  lengths must resolve to distinct bins.
- `with_length_range`: whole bins overlapping one half-open `[start, end)` bp
  range. Edge bins are not split or prorated.
- `length_bin_idxs`: direct bin indices.

Use only one length-bin selector at a time.

`value` means:

- `count`: loaded count value.
- `fraction`: count divided by the selected row's total count over length bins.
- `density`: `fraction` divided by `length_width_bp`, so unequal length-bin
  widths remain comparable.

`denominator` controls the row total used for `fraction` and `density`:

- `all_bins`: divide by the row total over all length bins.
- `selected_bins`: divide by the row total over the returned length bins.

`denominator` has no effect for `count`.

`max_blacklisted_fraction` is valid only for row modes that can carry
`blacklisted_fraction`. It must be a finite value in `[0, 1]`. The default
`1.0` keeps all rows when blacklist metadata is present. Stricter cutoffs error
clearly when the loaded output has no `blacklisted_fraction` column.

## Midpoint Profiles

Python:

```python
midpoints = cfl.read_midpoints("sample.midpoint_profiles.zarr")

midpoints.group_metadata()
midpoints.length_bins()
midpoints.positions()
midpoints.group_idx("LYL1")
midpoints.length_bin_idx(167)

midpoints.counts_array(
    groups=None,
    group_idxs=None,
    with_lengths=None,
    with_length_range=None,
    length_bin_idxs=None,
) # shape is always (group, length_bin, position)

midpoints.data_frame(
    groups=None,
    group_idxs=None,
    with_lengths=None,
    with_length_range=None,
    length_bin_idxs=None,
)
```

R:

```r
midpoints <- read_midpoints("sample.midpoint_profiles.zarr")

group_metadata(midpoints)
length_bins(midpoints)
positions(midpoints)
group_idx(midpoints, "LYL1")
length_bin_idx(midpoints, 167)

midpoint_array(midpoints)
profile_array(midpoints, group_idx = 1L, length_bin_idx = 1L)

midpoint_data_frame(
  midpoints,
  groups = NULL,
  group_idxs = NULL,
  with_lengths = NULL,
  with_length_range = NULL,
  length_bin_idxs = NULL
)
```

`with_lengths` selects the length bins containing the requested fragment
lengths. Multiple `with_lengths` values must resolve to distinct length bins;
if two lengths fall in the same bin, pass one representative length or use direct
bin indices. `with_length_range` selects whole bins overlapping a half-open
`[start, end)` bp range. `length_bin_idxs` selects bins directly. Use only one
length-bin selector at a time.

## End-Motif Counts

Python:

```python
ends = cfl.read_end_motifs("sample.end_motifs.zarr")

ends.storage_mode()
ends.motifs_metadata()
ends.motif_idx("_AA")
ends.has_motif("_AA")

ends.dense_counts_array(
    motifs=None,
    motif_idxs=None,
    allow_densify=False,
) # global, shape is always (1, motif)
ends.sparse_counts_matrix(
    motifs=None,
    motif_idxs=None,
) # global, shape is always (1, motif)
ends.window_metadata() # windowed only
ends.group_metadata()  # grouped only

ends.dense_counts_array(
    window_idxs=None,
    motifs=None,
    motif_idxs=None,
    allow_densify=False,
) # windowed, shape is always (window, motif)
ends.sparse_counts_matrix(
    window_idxs=None,
    motifs=None,
    motif_idxs=None,
) # windowed, shape is always (window, motif)

ends.dense_counts_array(
    groups=None,
    group_idxs=None,
    motifs=None,
    motif_idxs=None,
    allow_densify=False,
) # grouped, shape is always (group, motif)
ends.sparse_counts_matrix(
    groups=None,
    group_idxs=None,
    motifs=None,
    motif_idxs=None,
) # grouped, shape is always (group, motif)

ends.data_frame(
    densify=False,
    motifs=None,
    motif_idxs=None,
) # global
ends.data_frame(
    window_idxs=None,
    densify=False,
    motifs=None,
    motif_idxs=None,
    max_blacklisted_fraction=1.0,
) # windowed
ends.data_frame(
    groups=None,
    group_idxs=None,
    densify=False,
    motifs=None,
    motif_idxs=None,
    max_blacklisted_fraction=1.0,
) # grouped
```

R:

```r
ends <- read_end_motifs("sample.end_motifs.zarr")

storage_mode(ends)
motifs(ends)
motif_idx(ends, "_AA")
has_motif(ends, "_AA")

dense_counts_vector(ends) # global dense-compatible output
dense_counts_matrix(ends) # windowed or grouped dense-compatible output
sparse_counts_matrix(ends)

end_motif_data_frame(
  ends,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL
) # global
end_motif_data_frame(
  ends,
  window_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0
) # windowed
end_motif_data_frame(
  ends,
  groups = NULL,
  group_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0
) # grouped
```

For sparse stores, `densify = FALSE` returns only stored non-zero motif-count
rows. `densify = TRUE` adds explicit zero-count rows for the selected rows and
the selected observed motifs. It does not invent rows for motifs absent from the
file. Dense stores ignore `densify` and return their available dense grid.

Dense and sparse end-motif data frames use the same public columns for a given
row mode. Only the set of returned rows differs.

Schema version 2 end-motif stores may mark the motif axis as `motif_group`.
R and Python loaders still expose that count-column axis as `motifs` and
`motif_idxs`. In that case the returned motif labels are user-defined group
names, not concrete DNA motifs.
