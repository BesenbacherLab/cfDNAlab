# cfdnalab

Python helpers for loading cfDNAlab output files.

This package does not install the cfDNAlab command-line tool. The CLI is
distributed separately as the Rust `cfdna` binary. Use this Python package after
running cfDNAlab to load and analyze output files.

The first supported output type is midpoint Zarr output:
`<prefix>.midpoint_profiles.zarr`.

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

## Inspect Metadata

```python
groups = midpoints.groups()
length_bins = midpoints.length_bins()
positions = midpoints.positions()
```

`groups()` returns `group_idx`, `group_name`, and `eligible_intervals`.
`length_bins()` and `positions()` return the corresponding bin indices and
half-open bp coordinates.

## Extract One Profile

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

## Filter By Eligible Intervals

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

## NumPy Arrays

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

## Test

```bash
cd py_cfdnalab
uv run pytest
```
