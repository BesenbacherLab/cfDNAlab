# Improved Analysis Output Formats Plan

This plan covers a future cleanup of public analysis outputs. The goal is to
make outputs easier to trust and use from R and Python without forcing users to
reconstruct meaning from a raw NumPy array plus loosely related sidecars.

The main correction from the earlier Zarr-only idea is that not every output
needs Zarr. `lengths` is always two-dimensional and is naturally tabular.
`midpoints` and sparse/dense motif outputs are more array-like and are better
Zarr candidates.

## Format Direction

- `lengths`: write a wide, compressed TSV as the main user-facing output.
- `midpoints`: write a Zarr store with named dimensions and coordinate arrays.
- `ends`: write Zarr for dense and sparse motif outputs, including sparse COO
  arrays inside the Zarr store instead of SciPy-style `.npz`.
- Internal temporary files and GC correction packages can stay as they are until
  there is a separate reason to change them.

This keeps simple outputs simple and uses Zarr where it solves a real structure
problem.

## Lengths Output

`cfdna lengths` currently writes a dense matrix with shape `(row, length_bin)`.
That shape is always two-dimensional:

- One row for global mode.
- One row per genomic window for fixed-window and BED modes.
- One row per group for grouped-BED mode.
- One count column per length bin.

A long tidy TSV would be easy to plot, but it would repeat row metadata for every
length bin and can become very large. A wide TSV avoids that repetition while
remaining easy to load in R and Python.

### Wide TSV Contract

Write:

```text
<prefix>.length_counts.tsv.gz
<prefix>.fragment_length_settings.json
```

The TSV should use count columns named:

```text
count_<length>
count_<start>_<end>
```

Single-bp bins use `count_<length>` and mean:

```text
[length, length + 1)
```

Wider bins use `count_<start>_<end>` and mean:

```text
[start, end)
```

No separate `length_bins.tsv` is needed when the bin edges are encoded in the
count column names. Do not store midpoint base pairs. Users can derive plotting
midpoints if they need them.

### Global Mode

Global output has one row and no row key:

```text
count_30  count_31  count_32
123       140       131
```

### Fixed-Window And BED Modes

Non-grouped window output should be keyed only by coordinates:

```text
chrom  start  end  blacklisted_fraction  count_30  count_31  count_32
chr1   0      100000  0.001              123       140       131
chr1   100000 200000  0.000              55        61        70
```

Do not expose `row_index`. For BED mode, preserve current output ordering by the
original BED order, but do not require users to join on that index.

### Grouped BED Mode

Grouped output should be keyed only by `group_name`:

```text
group_name  eligible_windows  blacklisted_fraction  count_30  count_31  count_32
TSS         18342             0.012                 123       140       131
Enhancer    9481              0.009                 55        61        70
```

`eligible_windows` is the number of grouped BED windows retained in the group
after chromosome filtering and any command-level window filtering that affects
the grouped output. It is useful as a denominator when users want mean length
profiles per eligible window:

```text
mean_count_* = count_* / eligible_windows
```

This is not currently written by `lengths`, but the grouped-window map contains
enough information to compute it. Use the name `eligible_windows`, not
`eligible_intervals`, because these are length-count window assignment units.

`blacklisted_fraction` should keep its current grouped meaning: interval-width
weighted blacklisted fraction across windows in the group. If no blacklist was
used, either omit the column or write `0.0` consistently. Omitting it matches the
current grouped sidecar behavior more closely.

### Loading Lengths In R

```r
length_counts <- read.delim("sample.length_counts.tsv.gz")
```

For grouped output:

```r
count_columns <- grep("^count_[0-9]+(_[0-9]+)?$", names(length_counts), value = TRUE)

long_counts <- reshape(
  length_counts,
  varying = count_columns,
  v.names = "count",
  timevar = "length_bin",
  times = count_columns,
  direction = "long"
)
```

Users can parse `length_bin` names only when they need numeric starts and ends.

### Loading Lengths In Python

```python
import pandas as pd

length_counts = pd.read_csv("sample.length_counts.tsv.gz", sep="\t")
count_columns = [column for column in length_counts.columns if column.startswith("count_")]

long_counts = length_counts.melt(
    id_vars=[column for column in length_counts.columns if column not in count_columns],
    value_vars=count_columns,
    var_name="length_bin",
    value_name="count",
)

edges = long_counts["length_bin"].str.extract(r"count_(\d+)(?:_(\d+))?")
long_counts["length_start_bp"] = edges[0].astype(int)
long_counts["length_end_bp"] = edges[1].fillna(
    long_counts["length_start_bp"] + 1
).astype(int)
```

This keeps the file compact by default and still gives users a direct path to a
plotting dataframe.

## Midpoints Output

`cfdna midpoints` should be the first Zarr target. Its primary output is a real
three-dimensional tensor:

```text
counts[group, length_bin, position]
```

The current `.npy` plus `group_index.tsv` plus settings JSON makes users piece
together the group axis, length axis, and position axis manually. Zarr is useful
here because the dimensions and coordinates belong with the array.

### Midpoint Zarr Store

Write:

```text
<prefix>.midpoint_profiles.zarr/
```

Core arrays:

```text
counts:                  float32[group, length_bin, position]
group_name:              string[group]
eligible_intervals:      uint64[group]
length_start_bp:         uint32[length_bin]
length_end_bp:           uint32[length_bin]
position_start_bp:       int32[position]
position_end_bp:         int32[position]
```

Root attrs:

```json
{
  "cfdnalab_schema": "midpoint_profiles",
  "cfdnalab_schema_version": 1,
  "primary_array": "counts",
  "count_units": "weighted_midpoint_count",
  "settings": {
    "...": "existing midpoint_profile_settings.json content"
  }
}
```

Do not make users rely on implicit axis order. The array metadata should name
dimensions as:

```text
["group", "length_bin", "position"]
```

### Midpoints In Python

```python
import xarray as xr

dataset = xr.open_zarr("sample.midpoint_profiles.zarr", consolidated=False)

profile_table = (
    dataset["counts"]
    .isel(group=slice(0, 10))
    .to_dataframe(name="count")
    .reset_index()
)
```

Expected dataframe columns:

```text
group
length_bin
position
count
group_name
eligible_intervals
length_start_bp
length_end_bp
position_start_bp
position_end_bp
```

### Midpoints In R

```r
library(zarr)

store <- open_zarr("sample.midpoint_profiles.zarr", read_only = TRUE)

profiles <- store[["/counts"]][1:10, , ]
group_name <- store[["/group_name"]][]
length_start_bp <- store[["/length_start_bp"]][]
length_end_bp <- store[["/length_end_bp"]][]
position_start_bp <- store[["/position_start_bp"]][]
position_end_bp <- store[["/position_end_bp"]][]
```

For large 3D arrays, docs should explicitly tell users to subset before building
a dataframe. A full 3D dataframe can be much larger than the Zarr array.

## Ends Output

`cfdna ends` currently has two user-facing matrix modes:

- Sparse observed motifs as `.sparse.npz` plus motif labels.
- Dense all-motif output as `.npy` plus motif labels.

If Zarr is adopted for this command, do not keep a mixed `.npz` plus Zarr story.
Use Zarr for both dense and sparse outputs.

### Dense Ends Zarr

For dense all-motif output:

```text
<prefix>.end_motifs.zarr/
```

Arrays:

```text
counts:                  float64[row, motif]
motif:                   string[motif]
```

Row metadata follows the same command mode rules as other windowed outputs:

- Global mode has one row and no row key.
- Non-grouped window mode uses `chrom`, `start`, `end`.
- Grouped mode uses `group_name`.

Dimension names:

```text
["row", "motif"]
```

### Sparse Ends Zarr

For sparse output, store a COO representation inside the Zarr group:

```text
<prefix>.end_motifs.zarr/
```

Arrays:

```text
sparse/row:              uint64[nnz]
sparse/col:              uint64[nnz]
sparse/data:             float64[nnz]
sparse/shape:            uint64[2]
motif:                   string[motif]
```

Attrs:

```json
{
  "storage": "sparse_coo",
  "logical_shape": ["row", "motif"],
  "index_base": 0
}
```

This is enough to reconstruct sparse matrices:

Python:

```python
import scipy.sparse
import zarr

root = zarr.open_group("sample.end_motifs.zarr", mode="r")
matrix = scipy.sparse.coo_matrix(
    (root["sparse/data"][:], (root["sparse/row"][:], root["sparse/col"][:])),
    shape=tuple(root["sparse/shape"][:]),
)
motifs = root["motif"][:]
```

R:

```r
library(Matrix)
library(zarr)

store <- open_zarr("sample.end_motifs.zarr", read_only = TRUE)

row <- store[["/sparse/row"]][] + 1L
col <- store[["/sparse/col"]][] + 1L
data <- store[["/sparse/data"]][]
shape <- store[["/sparse/shape"]][]

matrix <- sparseMatrix(
  i = row,
  j = col,
  x = data,
  dims = as.integer(shape)
)
```

A sparse Zarr group is not a universal standard sparse format, but it is
transparent and reconstructable. It avoids keeping motif data in Zarr while
leaving sparse coordinates in NPZ.

## Shared Design Rules

- Do not expose row indices in public TSVs when domain keys exist.
- Use `group_name` as the grouped output key.
- Use genomic coordinates as the non-grouped window key.
- Encode half-open single-bp length bins as `count_<length>` in wide TSVs.
- Encode wider half-open length bins as `count_<start>_<end>` in wide TSVs.
- Do not store midpoint base pairs for length bins.
- Store coordinate arrays in Zarr for multidimensional outputs.
- Put settings and provenance in attrs or JSON sidecars, but do not make attrs
  the only place where axis coordinates live.
- Keep sparse arrays explicit: `row`, `col`, `data`, and `shape`.
- Use zero-based sparse coordinates on disk. R examples should add 1 when
  constructing R sparse matrices.

## Implementation Order

1. Add wide TSV output for `lengths`. Implemented in the current working tree.
2. Add `eligible_windows` for grouped `lengths`. Implemented in the current working tree.
3. Add Zarr output for `midpoints`.
4. Move `ends` dense output to Zarr.
5. Move `ends` sparse output to sparse COO Zarr.
6. Update website docs with R and Python loading examples.
7. Only then decide whether positional k-mer and transition outputs should join
   the same Zarr pattern.

## Open Decisions

- Whether `blacklisted_fraction` should be omitted or written as `0.0` when no
  blacklist was used.
- Whether Zarr V3 should be required immediately, or whether V2 is needed for R
  compatibility.
- Which codec is safest across Python and R readers.
- Whether sparse Zarr should use COO only, or also include CSR-style `indptr` for
  faster row access.
- Whether the output directory should contain both TSV and Zarr outputs for some
  commands, or whether each command should have one primary machine-readable
  format.

## Acceptance Criteria

The transition is successful when:

- `lengths` can be loaded directly with `read.delim()` or `pandas.read_csv()`.
- `lengths` grouped output exposes `group_name` and `eligible_windows`.
- `lengths` non-grouped output exposes genomic coordinates without row indices.
- `midpoints` opens as a named, coordinate-aware dataset in Python and R.
- `ends` sparse output can be reconstructed from Zarr in Python and R.
- Users no longer have to manually combine a raw array, motif labels, group
  labels, and settings just to make a trustworthy analysis dataframe.
