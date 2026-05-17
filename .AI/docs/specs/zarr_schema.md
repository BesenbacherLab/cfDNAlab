# Zarr Schema Spec

cfDNAlab public analysis arrays are Zarr V3 stores. Zarr metadata owns chunking,
compression, dtype encoding, and native `dimension_names`; cfDNAlab root
attributes and array attributes define the biological schema.

## Shared Rules

- Store paths must end in `.zarr` and contain a root `zarr.json`.
- Root attributes must include `cfdnalab_schema` and
  `cfdnalab_schema_version`.
- Current schema version is `1` for all Zarr stores in this spec.
- `int32` arrays are used for coordinate axes, small label indices, and small
  non-negative metadata that should load as native R integers.
- `uint64` arrays are used for genomic coordinates.
- R readers should preserve `uint64` as `bit64::integer64` or check that values
  are below `2^53` before conversion to double. Sparse indices are `int32` and
  must be range checked before conversion to R's one-based matrix indices.
- Public arrays use fill values outside their valid data domain where possible.
  This prevents readers that map Zarr fill values to missing values from
  turning valid zero counts or zero-based coordinates into missing values.
  Current fill values are `-1` for non-negative `int32` metadata, `-1.0` for
  non-negative floating counts/fractions, `u64::MAX` for genomic coordinates,
  and `255` for fixed-width ASCII label bytes. The `255` fill value must not be
  reused for arbitrary numeric `uint8` arrays where `255` could be real data.

## Midpoint Profiles

Root attributes:

- `cfdnalab_schema = "midpoint_profiles"`
- `cfdnalab_schema_version = 1`
- `primary_array = "counts"`
- `count_units = "weighted_midpoint_count"`

Arrays:

- `counts[group, length_bin, position]`: `float32`, primary weighted midpoint
  count tensor.
- `group[group]`: `int32`, zero-based group index. Attributes:
  `label_field = "group_name"` and `labels = [...]`.
- `eligible_intervals[group]`: `int32`, retained profile intervals per group.
- `length_bin[length_bin]`: `int32`, zero-based length-bin index.
- `length_start_bp[length_bin]`: `int32`, inclusive fragment length-bin start.
- `length_end_bp[length_bin]`: `int32`, exclusive fragment length-bin end.
- `position[position]`: `int32`, zero-based profile position-bin index.
- `position_bin_start_bp[position]`: `int32`, inclusive interval-relative
  position-bin start.
- `position_bin_end_bp[position]`: `int32`, exclusive interval-relative
  position-bin end.

The root does not duplicate `dimension_names`. Readers should use the native
Zarr V3 array metadata.

## End-Motif Counts

Root attributes:

- `cfdnalab_schema = "end_motif_counts"`
- `cfdnalab_schema_version = 1`
- `storage_mode = "dense"` or `"sparse_coo"`
- `row_mode = "global"`, `"size"`, `"bed"`, or `"grouped_bed"`
- `count_units = "weighted_end_motif_count"`
- `primary_array = "counts"` and `primary_group = null` for dense output.
- `primary_array = null`, `primary_group = "sparse"`,
  `sparse_format = "coo"`, and `sparse_indices_base = 0` for sparse output.

Motif metadata:

- `motif_index[motif]`: `int32`, zero-based motif column index.
- `motif_byte[motif_byte]`: `int32`, zero-based byte offset within each motif
  label.
- `motif_ascii[motif, motif_byte]`: `uint8`, fixed-width ASCII motif labels.
  Decode each motif row as ASCII to recover labels in count-column order.

Shared row metadata:

- `row[row]`: `int32`, zero-based count-row index. Global output adds
  `label_field = "row_label"` and `labels = ["global"]`.

Windowed row metadata for `row_mode = "size"` or `"bed"`:

- `chromosome[chromosome]`: `int32`, chromosome index. Attributes:
  `label_field = "chromosome_name"` and `labels = [...]`.
- `row_chromosome[row]`: `int32`, index into `chromosome`.
- `row_start_bp[row]`: `uint64`, inclusive row start coordinate.
- `row_end_bp[row]`: `uint64`, exclusive row end coordinate.
- `blacklisted_fraction[row]`: `float64`, row blacklist fraction.

Grouped row metadata for `row_mode = "grouped_bed"`:

- `group[row]`: `int32`, group index matching the count row. Attributes:
  `label_field = "group_name"` and `labels = [...]`.
- `eligible_windows[row]`: `int32`, retained grouped BED windows per group.
- `blacklisted_fraction[row]`: `float64`, length-weighted group blacklist
  fraction.

Dense counts:

- `counts[row, motif]`: `float64`, primary weighted end-motif count matrix.

Sparse counts:

- `sparse/row[nnz]`: `int32`, zero-based COO row indices.
- `sparse/motif[nnz]`: `int32`, zero-based COO motif indices.
- `sparse/count[nnz]`: `float64`, weighted non-zero counts.
- `sparse/shape[sparse_dimension]`: `int32`, dense shape `[row, motif]`.
- `sparse/sparse_dimension[sparse_dimension]`: `int32`, axis index with
  `label_field = "sparse_dimension_name"` and `labels = ["row", "motif"]`.

Sparse COO entries are sorted by `(row, motif)` and must not contain duplicate
coordinate pairs.
