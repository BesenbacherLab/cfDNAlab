# Zarr Schema Spec

cfDNAlab public analysis arrays are Zarr V3 stores. Zarr metadata owns chunking,
compression, dtype encoding, and native `dimension_names`; cfDNAlab root
attributes and array attributes define the biological schema.

## Shared Rules

- Store paths must end in `.zarr` and contain a root `zarr.json`.
- Root attributes must include `cfdnalab_schema` and
  `cfdnalab_schema_version`.
- Current schema versions are `1` for `midpoint_profiles` and
  `end_motif_counts`, and `3` for `reference_gc_package` and
  `gc_correction_package`.
- `int32` arrays are used for coordinate axes, small label indices, and small
  non-negative metadata that should load as native R integers.
- `int64` arrays are used for genomic coordinates. They keep the coordinate
  domain large while avoiding `uint64` handling problems in R Zarr readers.
  Sparse indices are `int32` and must be range checked before conversion to R's
  one-based matrix indices.
- Public arrays use fill values outside their valid data domain where possible.
  This prevents readers that map Zarr fill values to missing values from
  turning valid zero counts or zero-based coordinates into missing values. The
  midpoint and end-motif schemas use `-1` for non-negative `int32` metadata,
  `-1.0` for non-negative floating counts/fractions, `-1` for genomic
  coordinates, and `255` for fixed-width ASCII label bytes. The `255` fill
  value must not be reused for arbitrary numeric `uint8` arrays where `255`
  could be real data. GC package arrays are currently single-chunk package
  arrays and use domain-specific zero or false fill values.

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
- `row_start_bp[row]`: `int64`, inclusive row start coordinate.
- `row_end_bp[row]`: `int64`, exclusive row end coordinate.
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

## Reference GC Package

Root attributes:

- `cfdnalab_schema = "reference_gc_package"`
- `cfdnalab_schema_version = 3`
- `package_role = "reference_gc"`
- `value_units = "reference_fragment_mass"`
- `gc_percent_rounding = "integer_half_up"`
- `minimum_acgt_bases_for_gc_fraction`: minimum ACGT bases required before GC
  fraction is defined.
- `end_offset`, `skip_interpolation`, `smoothing_radius`, `smoothing_sigma`,
  and `skip_smoothing`: settings needed to validate and reuse the package.

Arrays:

- `counts[length, gc_percent]`: `float64`, expected reference fragment mass.
- `support_mask_unobservables[length, gc_percent]`: `bool`, theoretical
  reference support mask.
- `support_mask_outliers[length, gc_percent]`: `bool`, empirical reference
  support mask.
- `gc_percent_widths[length, gc_percent]`: `uint16`, number of integer GC
  counts represented by each GC-percent bin.
- `length[length]`: `int32`, fragment length in bp.
- `gc_percent[gc_percent]`: `int32`, integer GC-percent axis.
- `chromosome[chromosome]`: `int32`, selected chromosome index. Attributes:
  `label_field = "chromosome_name"` and `labels = [...]`.
- `reference_contig_footprint_json[json_byte]`: `uint8`, JSON-encoded
  reference contig footprint used for compatibility checks.

## GC Correction Package

Root attributes:

- `cfdnalab_schema = "gc_correction_package"`
- `cfdnalab_schema_version = 3`
- `package_role = "sample_gc_correction"`
- `correction_units = "multiplicative_fragment_weight"`
- `gc_percent_rounding = "integer_half_up"`
- `minimum_acgt_bases_for_gc_fraction`: minimum ACGT bases required before GC
  fraction is defined.
- `end_offset`: fragment-end trimming offset used when computing GC.

Arrays:

- `correction_matrix[length_bin, gc_bin]`: `float64`, multiplicative GC
  correction weights.
- `length_edges[length_edge]`: `uint32`, half-open fragment length-bin edges.
- `gc_edges[gc_edge]`: `uint32`, half-open GC-percent bin edges.
- `length_bin_frequencies[length_bin]`: `float64`, normalized sample
  length-bin frequency weights.
- `reference_contig_footprint_json[json_byte]`: `uint8`, JSON-encoded
  reference contig footprint inherited from the reference GC package.
