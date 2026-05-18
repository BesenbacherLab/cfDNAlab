# cfDNAlab | Python Loaders <img src="https://raw.githubusercontent.com/BesenbacherLab/cfDNAlab/refs/heads/main/cfdnalab_logo_little_guy_172x200_144dpi.png" align="right" height="155" />

Python helpers for loading [**cfDNAlab**](https://github.com/BesenbacherLab/cfDNAlab) output files.

This package does not install or run the cfDNAlab command-line tool. The CLI is distributed separately as the Rust `cfdna` binary. Use this Python package after running cfDNAlab to load and analyze output files.

The first supported output types are midpoint and end-motif Zarr outputs: `<prefix>.midpoint_profiles.zarr` and `<prefix>.end_motifs.zarr`.

<br>

## Install

These instructions only installs the Python loader package. To install the `cfdna` command-line tool, see the [main repository](https://github.com/BesenbacherLab/cfDNAlab).

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

profile = midpoints.data_frame_for_profile(
    group_idx=group_idx,
    length_bin_idx=length_bin_idx,
)
```

The returned data frame has one row per midpoint position bin.

### Filter By Eligible Intervals

```python
min_intervals = 100

for _, group in midpoints.groups().iterrows():
    if group["eligible_intervals"] < min_intervals:
        continue

    profile = midpoints.data_frame_for_profile(
        group_idx=group["group_idx"],
        length_bin_idx=0,
    )
```

### Extract NumPy Arrays

```python
profile = midpoints.array_for_profile(group_idx=0, length_bin_idx=0)
group_counts = midpoints.array_from_group_idx(group_idx=0)
length_counts = midpoints.array_from_length_bin(length_bin_idx=0)
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

For sparse output, `sparse_coo_data_frame()` is usually the easiest way to inspect or plot the non-zero motif counts. Use `sparse_coo()` or the sparse slice helpers when you want SciPy sparse matrices. Dense helpers require `allow_densify=True` on sparse stores so large sparse outputs are not accidentally expanded in memory.

For dense output, the `dense_data_frame*()` methods are usually the most convenient starting point. Use `dense_counts_zarr_array()` when you want the on-disk Zarr array and `dense_counts_matrix()` when you want the full NumPy matrix in memory.

`sparse_coo_data_frame()` is only available for sparse output.

### Inspect End-Motif Metadata

```python
motifs = ends.motif_metadata()
ends.has_motif("_AA")
```

`read_end_motifs()` returns a mode-specific object.

- Windowed output has `windows()`, which returns `window_idx`, `chrom`, `start`,
  `end`, and `blacklisted_fraction`.
- Grouped output has `groups()` and `group_idx()`.
- Global output has `dense_counts_vec()` and `dense_data_frame()`.

### Extract End-Motif Counts

```python
motif_idx = ends.motif_idx("_AA")

motif_counts = ends.dense_data_frame_for_motif_idx(motif_idx)
```

Sparse output stays sparse unless you ask for dense arrays:

```python
sparse_counts = ends.sparse_coo()
sparse_payload = ends.sparse_coo_data_frame()
motif_array = ends.dense_counts_for_motif("_AA", allow_densify=True)
```

For dense windowed output:

```python
windows = ends.windows()
window_counts = ends.dense_data_frame_for_window(window_idx=0)
```

For dense grouped output:

```python
groups = ends.groups()
group_counts = ends.dense_data_frame_for_group("t-cells")
```

For sparse stores, prefer `sparse_coo()`, `sparse_coo_data_frame()`, and the sparse slice helpers when working with large end-motif outputs. Use `allow_densify=True` only when the dense result is small enough to fit comfortably in memory.
