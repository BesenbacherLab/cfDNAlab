# cfDNAlab <img src='https://raw.githubusercontent.com/BesenbacherLab/cfDNAlab/refs/heads/main/cfdnalab_logo_little_guy_172x200_144dpi.png' align="right" height="140" />

Python helpers for loading [**cfDNAlab**](https://github.com/BesenbacherLab/cfDNAlab) output files.

**NOTE**: This package does not install the cfDNAlab command-line tool for extracting features. 
The CLI is distributed separately as the Rust `cfdna` binary. Use this Python package after
running cfDNAlab to load and analyze output files.

The first supported output types are midpoint and end-motif Zarr outputs:
`<prefix>.midpoint_profiles.zarr` and `<prefix>.end_motifs.zarr`.

## Install

From this repository:

```bash
cd py_cfdnalab
uv sync --extra test
```

For an editable install in another environment:

```bash
cd py_cfdnalab
uv pip install -e .
```

## Load Midpoint Profiles

```python
import cfdnalab as cfl

midpoints = cfl.load_midpoints("sample.midpoint_profiles.zarr")
```

### Inspect Metadata

```python
groups = midpoints.groups()
length_bins = midpoints.length_bins()
positions = midpoints.positions()
```

`groups()` returns `group_idx`, `group_name`, and `eligible_intervals`.
`length_bins()` and `positions()` return the corresponding bin indices and
half-open bp coordinates.

### Extract One Profile

Use `group_idx()` and `length_bin()` when selecting by names or bp lengths:

```python
group_idx = midpoints.group_idx("CTCF")
length_bin = midpoints.length_bin(167)

profile = midpoints.data_frame_for_profile(
    group_idx=group_idx,
    length_bin=length_bin,
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
        length_bin=0,
    )
```

### Extract NumPy Arrays

```python
profile = midpoints.array_for_profile(group_idx=0, length_bin=0)
group_counts = midpoints.array_from_group_idx(group_idx=0)
length_counts = midpoints.array_from_length_bin(length_bin=0)
```

`array()` loads the full 3D count tensor into RAM:

```python
counts = midpoints.array()
```

Prefer the slice helpers when possible.

## Load End-Motif Counts

```python
import cfdnalab as cfl

ends = cfl.load_end_motifs("sample.end_motifs.zarr")
```

### Storage Mode - Sparse or Dense

Start by checking whether the counts were stored as a dense matrix or sparse COO arrays.

```python
ends.storage_mode()
```

If the storage mode is `"sparse_coo"`, the `sparse_coo*()` methods use the stored COO arrays without densifying. Use the `dense_*()` methods when you explicitly want dense NumPy arrays or dense data frames.

If the storage mode is `"dense"`, `sparse_coo*()` methods still work, but they first read the dense count matrix and convert it to a SciPy sparse object.

### Inspect End-Motif Metadata

```python
motifs = ends.motif_metadata()
```

`load_end_motifs()` returns a mode-specific object. Windowed output has `windows()`. Grouped output has `groups()` and `group_idx()`. Global output has `counts()` and `data_frame()`.

### Extract End-Motif Counts

```python
motif_idx = ends.motif_idx("_AA")

motif_counts = ends.dense_data_frame_for_motif_idx(motif_idx)
```

Sparse output stays sparse unless you ask for dense arrays:

```python
sparse_counts = ends.sparse_coo()
sparse_payload = ends.sparse_coo_data_frame()
motif_array = ends.dense_array_for_motif("_AA")
```

For windowed output:

```python
windows = ends.windows()
window_counts = ends.dense_data_frame_for_window(window_idx=0)
```

For grouped output:

```python
groups = ends.groups()
group_counts = ends.dense_data_frame_for_group("CTCF")
```

Methods prefixed with `dense_` may densify sparse output. Prefer `sparse_coo()`, `sparse_coo_data_frame()`, and the sparse slice helpers when working with large end-motif outputs.

## Test

```bash
cd py_cfdnalab
uv run pytest
```
