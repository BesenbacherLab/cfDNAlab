# cfDNAlab | Python Loaders <img src="https://raw.githubusercontent.com/BesenbacherLab/cfDNAlab/refs/heads/main/cfdnalab_logo_little_guy_172x200_144dpi.png" align="right" height="155" />

Python helpers for loading [**cfDNAlab**](https://github.com/BesenbacherLab/cfDNAlab) output files.

This package does not install or run the cfDNAlab command-line tool. The CLI is distributed separately as the Rust `cfdna` binary. Use this Python package after running cfDNAlab to load and analyze output files.

Supported output types are midpoint and end-motif Zarr outputs plus length-count TSV outputs: `<prefix>.midpoint_profiles.zarr`, `<prefix>.end_motifs.zarr`, and `<prefix>.length_counts.tsv.zst`.

<br>

## Install

These instructions install only the Python loader package. To install the `cfdna` command-line tool, see the [main repository](https://github.com/BesenbacherLab/cfDNAlab).

Install with pip:

```bash
pip install cfdnalab
```

Install the current development version from GitHub:

```bash
pip install "cfdnalab @ git+https://github.com/BesenbacherLab/cfDNAlab.git#subdirectory=py-cfdnalab"
```

<br>

## Load Midpoint Profiles

```python
import cfdnalab as cfl

midpoints = cfl.read_midpoints("sample.midpoint_profiles.zarr")
```

### Inspect Metadata

```python
groups = midpoints.groups()
length_bins = midpoints.length_bins()
positions = midpoints.positions()
```

`groups()` returns `group_idx`, `group_name`, and `eligible_intervals`. `length_bins()` and `positions()` return the corresponding bin indices and half-open bp coordinates.

### Extract One Profile

Use `group_idx()` and `length_bin_idx()` when selecting by names or bp lengths:

```python
group_idx = midpoints.group_idx("LYL1")
length_bin_idx = midpoints.length_bin_idx(167)

profile = midpoints.data_frame(group_idxs=group_idx, length_bin_idxs=length_bin_idx)
```

The returned data frame has one row per midpoint position bin.

### Extract A Group Or Length Bin

Use `data_frame(groups=...)` for all length and position bins in one group. Use `data_frame(with_lengths=...)` when you have a fragment length in bp and want the length bin that contains it.

```python
group_data = midpoints.data_frame(groups="LYL1")
length_bin_data = midpoints.data_frame(with_lengths=167)
```

### Filter By Eligible Intervals

```python
min_intervals = 100

for _, group in midpoints.groups().iterrows():
    if group["eligible_intervals"] < min_intervals:
        continue

    profile = midpoints.data_frame(group_idxs=group["group_idx"], length_bin_idxs=0)
```

### Extract NumPy Arrays

```python
profile = midpoints.array_for_profile(group_idx=0, length_bin_idx=0)
group_counts = midpoints.array_from_group("LYL1")
length_bin_counts = midpoints.array_from_length(167)
```

`array()` loads the full 3D count tensor into RAM:

```python
counts = midpoints.array()
```

Prefer the slice helpers when possible.

<br>

## Load End-Motif Counts

```python
import cfdnalab as cfl

ends = cfl.read_end_motifs("sample.end_motifs.zarr")
```

### Storage Mode - Sparse or Dense

Start by checking whether the counts were stored as a dense matrix or sparse COO arrays.

```python
ends.storage_mode()
```

For sparse output, `data_frame()` returns stored non-zero motif counts by default. Use `sparse_coo()` or the sparse slice helpers when you want SciPy sparse matrices. Pass `densify=True` only when the zero-filled result is small enough to fit in memory. Densifying only includes observed motifs.

For dense output, `data_frame()` returns all selected rows and motifs. Use `dense_counts_zarr_array()` when you want the on-disk Zarr array and `dense_counts_matrix()` when you want the full NumPy matrix in memory.

### Inspect End-Motif Metadata

```python
motifs = ends.motif_metadata()
motif_idx = ends.motif_idx("_AA")
ends.has_motif("_AA")
```

`read_end_motifs()` returns a mode-specific object.

- Windowed output has `windows()`, which returns `window_idx`, `chrom`, `start`,
  `end`, and `blacklisted_fraction`.
- Grouped output has `groups()` and `group_idx()`.
- Global output has `dense_counts_vec()` and `data_frame()`.

### Extract End-Motif Counts

```python
motif_idx = ends.motif_idx("_AA")

motif_counts = ends.data_frame(motifs="_AA")
```

Sparse output stays sparse unless you ask for dense arrays:

```python
nonzero_counts = ends.data_frame()
motif_count_matrix = ends.sparse_coo_for_motif("_AA")
motif_count_vector = ends.dense_counts_for_motif("_AA", allow_densify=True)
```

For dense windowed output:

```python
windows = ends.windows()
window_counts = ends.data_frame(window_idxs=0)
```

For dense grouped output:

```python
groups = ends.groups()
group_idx = ends.group_idx("t-cells")
group_counts = ends.data_frame(groups="t-cells")
```

For global output:

```python
global_counts = ends.dense_counts_vec(allow_densify=True)
global_data = ends.data_frame(densify=True)
```

For windowed or grouped outputs, `max_blacklisted_fraction` keeps rows with `blacklisted_fraction` at or below the cutoff:

```python
filtered_motif_counts = ends.data_frame(
    motifs="_AA",
    max_blacklisted_fraction=0.1,
)
```

For sparse stores, prefer `data_frame(densify=False)`, `sparse_coo()`, and the sparse slice helpers when working with large end-motif outputs. Use `densify=True` only when the dense result is small enough to fit comfortably in memory.

<br>

## Load Length Counts

```python
import cfdnalab as cfl

lengths = cfl.read_lengths("sample.length_counts.tsv.zst")
```

`read_lengths()` returns a mode-specific object. Windowed output has `windows()` and `counts_for_window()`, grouped output has `groups()`, `group_idx()`, and `counts_for_group()`, and global output has `counts_vec()`.

```python
bins = lengths.length_bins()
length_bin_idx = lengths.length_bin_idx(167)
counts = lengths.counts_matrix()
count_data = lengths.data_frame(value="count")
fraction_data = lengths.data_frame(value="fraction")
density_data = lengths.data_frame(value="density")
wide_density_data = lengths.data_frame(value="density", keep_wide=True)
```

For global output:

```python
global_counts = lengths.counts_vec()
global_data = lengths.data_frame(value="fraction")
```

For windowed output:

```python
windows = lengths.windows()
window_counts = lengths.counts_for_window(window_idx=0)
window_data = lengths.data_frame(window_idxs=0, value="fraction")
selected_windows = lengths.data_frame(
    window_idxs=[0, 2, 3],
    value="density",
    keep_wide=True,
)
```

For windowed or grouped outputs, `max_blacklisted_fraction` filters selected output rows before counts are returned:

```python
filtered = lengths.data_frame(max_blacklisted_fraction=0.1)
```

Outputs without a `blacklisted_fraction` column keep all rows at the default `max_blacklisted_fraction=1.0`. Stricter cutoffs raise an error as there is no blacklist column to filter on.

For grouped output:

```python
groups = lengths.groups()
group_idx = lengths.group_idx("t-cells")
group_counts = lengths.counts_for_group("t-cells")
group_data = lengths.data_frame(groups="t-cells", value="fraction")
selected_groups = lengths.data_frame(
    groups=["t-cells", "b-cells"],
    value="density",
    keep_wide=True,
    max_blacklisted_fraction=0.1,
)
```
