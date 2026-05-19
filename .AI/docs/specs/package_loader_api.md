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
lengths.counts_vector()      # global only
lengths.counts_matrix()      # windowed or grouped

lengths.data_frame(value="count", keep_wide=False)  # global
lengths.data_frame(
    window_idxs=None,
    value="count",
    keep_wide=False,
    max_blacklisted_fraction=1.0,
)  # windowed
lengths.data_frame(
    groups=None,
    group_idxs=None,
    value="count",
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
length_counts_matrix(lengths) # windowed or grouped

length_data_frame(lengths, value = "count", keep_wide = FALSE) # global
length_data_frame(
  lengths,
  window_idxs = NULL,
  value = "count",
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0
) # windowed
length_data_frame(
  lengths,
  groups = NULL,
  group_idxs = NULL,
  value = "count",
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0
) # grouped
```

Length data frames use the selected output rows and length bins. Long output has
one row per selected output unit and length bin. Wide output has one row per
selected output unit and one value column per length bin.

`value` means:

- `count`: loaded count value.
- `fraction`: count divided by the selected row's total count over length bins.
- `density`: `fraction` divided by `length_width_bp`, so unequal length-bin
  widths remain comparable.

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

midpoints.array()
midpoints.array_for_profile(group_idx=0, length_bin_idx=0)
midpoints.array_from_group("LYL1")
midpoints.array_from_group_idx(0)
midpoints.array_from_length(167)
midpoints.array_from_length_bin(0)

midpoints.data_frame(
    groups=None,
    group_idxs=None,
    with_lengths=None,
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
  length_bin_idxs = NULL
)
```

`with_lengths` selects the length bins containing the requested fragment
lengths. The same containing-bin rule applies to Python `array_from_length()`.
`length_bin_idxs` selects bins directly. Use one of those selectors, not both.

## End-Motif Counts

Python:

```python
ends = cfl.read_end_motifs("sample.end_motifs.zarr")

ends.storage_mode()
ends.motifs()
ends.motif_idx("_AA")
ends.has_motif("_AA")

ends.dense_counts_vector() # global dense-compatible output
ends.dense_counts_matrix() # windowed or grouped dense-compatible output
ends.sparse_counts_matrix()

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
