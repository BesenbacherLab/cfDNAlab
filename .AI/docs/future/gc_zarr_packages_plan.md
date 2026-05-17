# GC Zarr Packages Plan

## Scope

This plan is only about public command outputs that are also public command
inputs:

- `cfdna ref-gc-bias`: reusable reference GC package
- `cfdna gc-bias`: sample-specific GC correction package

Do not convert internal temporary `.npz` files in this work. Do not convert
`gc-bias --save-intermediates` staged `.npy` files in this work. Those files are
debug and QC artifacts, not the main package contract users pass between
commands.

The reason to do this now is that GC packages sit in the middle of real user
pipelines. If the public package format changes later, users either rerun
`ref-gc-bias` and `gc-bias` or need conversion tooling. That is a bigger
disruption than changing one terminal analysis output.

## Previous Public Surface

Reference package:

```text
<prefix>.ref_gc_package.npz
```

Current arrays:

```text
counts
support_mask_unobservables
support_mask_outliers
gc_percent_widths
version
length_range
end_offset
skip_interpolation
smoothing_radius
smoothing_sigma
skip_smoothing
chromosomes_json
reference_contig_footprint_json
```

Sample correction package:

```text
<prefix>.gc_bias_correction.npz
```

Current arrays:

```text
correction_matrix
length_edges
gc_edges
version
end_offset
length_bin_frequencies
reference_contig_footprint_json
```

The Zarr transition is a replacement, not a dual-format compatibility layer.
Downstream commands should validate public GC packages as `.zarr` directories.

## Target Public Surface

Reference package:

```text
<prefix>.ref_gc_package.zarr
```

Sample correction package:

```text
<prefix>.gc_bias_correction.zarr
```

The package directories are Zarr V3 stores written through `zarrs` with the same
shared Zarr helper layer used by midpoint and ends outputs. Use zstd compression.

The public package format is a replacement for the `.npz` package contract, not
a second parallel output.

## Root Metadata

Reference root attributes:

```json
{
  "cfdnalab_schema": "reference_gc_package",
  "cfdnalab_schema_version": 1,
  "package_role": "reference_gc",
  "value_units": "reference_fragment_mass",
  "gc_percent_rounding": "integer_half_up",
  "minimum_acgt_bases_for_gc_fraction": 10
}
```

Sample root attributes:

```json
{
  "cfdnalab_schema": "gc_correction_package",
  "cfdnalab_schema_version": 1,
  "package_role": "sample_gc_correction",
  "correction_units": "multiplicative_fragment_weight",
  "gc_percent_rounding": "integer_half_up",
  "minimum_acgt_bases_for_gc_fraction": 10
}
```

Keep command settings JSON outside the Zarr package if needed for human-readable
run provenance. The package root attributes should be the machine-readable
schema contract and small discovery metadata.

## Reference Package Schema

Dimensions:

```text
length       one row per integer fragment length
gc_percent   integer GC percent, usually 0..100
chromosome   selected chromosome names
json_byte    byte offset for compact JSON metadata payloads
```

Arrays:

```text
counts                       float64[length, gc_percent]
support_mask_unobservables   bool[length, gc_percent]
support_mask_outliers        bool[length, gc_percent]
gc_percent_widths            uint16[length, gc_percent]
length                       int32[length]
gc_percent                   int32[gc_percent]
chromosome                   int32[chromosome]
reference_contig_footprint_json uint8[json_byte]
```

`length` stores the actual integer fragment lengths from `min_fragment_length`
through `max_fragment_length`, not just row indices. This removes the need for a
separate `length_range` array.

`gc_percent` stores the actual integer percent labels. It should be `0..100` for
the current flow.

`chromosome` stores a numeric axis with JSON attrs:

```json
{
  "label_field": "chromosome_name",
  "labels": ["chr1", "chr2"]
}
```

Scalar or small settings should be root attributes, not one-element arrays:

```text
end_offset
skip_interpolation
smoothing_radius
smoothing_sigma
skip_smoothing
```

`reference_contig_footprint_json` remains a byte array because the footprint is
structured metadata. It is not a coordinate axis and should not be flattened into
many ad hoc arrays unless there is a concrete downstream reason.

## Sample Correction Package Schema

Dimensions:

```text
length_bin
length_edge
gc_bin
gc_edge
json_byte
```

Arrays:

```text
correction_matrix              float64[length_bin, gc_bin]
length_edges                   uint32[length_edge]
gc_edges                       uint32[gc_edge]
length_bin_frequencies         float64[length_bin]
reference_contig_footprint_json uint8[json_byte]
```

Root attributes:

```text
end_offset
```

The existing readback rule for edges should be preserved: edges are
inclusive/exclusive, with the final edge treated as inclusive on readback.

Keep `length_edges` and `gc_edges` as arrays because they are axis definitions,
not settings. Their lengths must be `n_length_bins + 1` and `n_gc_bins + 1`.

## Loader Changes

Reference GC loader:

- Replace `read_reference_gc_package(path)` with a Zarr reader that validates
  `cfdnalab_schema == "reference_gc_package"`.
- Validate schema version against a package-specific constant.
- Validate all matrix shapes against the coordinate axes.
- Build `ReferenceGCMetadata` from root attrs, the `length` axis, chromosome
  labels, and `reference_contig_footprint_json`.
- Keep all existing semantic validation after decoding:
  - non-empty chromosome labels
  - ordered length axis
  - effective minimum GC length >= 10
  - support masks and width table shape match counts
  - smoothing settings are valid when smoothing is enabled

Sample correction loader:

- Replace `GCCorrectionPackage::from_file` internals with a Zarr reader.
- Update path validation to require an existing `.zarr` directory.
- Validate `cfdnalab_schema == "gc_correction_package"`.
- Validate schema version against a package-specific constant.
- Preserve current correction-package validation:
  - edge lengths match correction matrix shape
  - length-bin frequencies match correction rows
  - correction weights are finite and non-negative
  - reference contig footprint is present and valid

Downstream `--gc-file` help and validation should say `.zarr`, not `.npz`.

## Writer Changes

Reference writer:

- Add a dedicated `src/commands/ref_gc_bias/zarr.rs` writer or move shared GC
  package writing into `src/commands/gc_bias/package_zarr.rs` if that keeps
  ownership clearer.
- Write the package store through `FinalOutputFiles`.
- Replace final path construction with:

```text
dot_join(&[prefix, "ref_gc_package.zarr"])
```

Sample writer:

- Replace `write_npz` with `GCCorrectionPackage::write_zarr`.
- Replace final path construction with:

```text
dot_join(&[prefix, "gc_bias_correction.zarr"])
```

Do not change cross-tile partial `.npz` files as part of this work.

Do not change `--save-intermediates` `.npy` files as part of this work.

## Reuse

Use shared Zarr helpers for:

- store creation
- root metadata writing
- dense array creation
- single-chunk coordinate arrays
- label validation for chromosome names
- final-output directory replacement safeguards already used for Zarr stores

Do not introduce generic command-agnostic schema builders. The GC package
schemas are small enough that explicit writers will be easier to review.

Useful small shared helpers may be added if they stay low-level:

- JSON byte-array writer
- JSON byte-array reader
- numeric coordinate-axis validation

## Python And R Downstream Helpers

Python helper package should eventually expose:

```python
ref_gc = cfdnalab.load_reference_gc_package("hg38.ref_gc_package.zarr")
gc = cfdnalab.load_gc_correction_package("sample.gc_bias_correction.zarr")
```

Minimum Python helper methods:

```text
ReferenceGCPackage:
  counts()
  support_masks()
  gc_percent_widths()
  lengths()
  gc_percents()
  chromosomes()

GCCorrectionPackage:
  correction_matrix()
  length_edges()
  gc_edges()
  length_bin_frequencies()
  data_frame()
```

R helpers should mirror those concepts with ordinary data frames and matrices.
The helpers are useful, but the core package must still be readable directly
with ordinary Zarr readers.

## Tests

Rust tests:

- Reference package writer round-trips through the new Zarr loader.
- Sample correction package writer round-trips through the new Zarr loader.
- Shape mismatches fail with clear errors.
- Missing root schema or wrong schema version fails with clear errors.
- Path validation rejects files and non-`.zarr` directories.
- Existing downstream correction tests continue to exercise
  `GCCorrectionPackage::from_file` through actual `.zarr` packages.

Downstream tests:

- Python `zarr` can read both stores directly.
- Python helper can load both stores and produce expected data frames/matrices.
- R Zarr reader can read both stores directly.
- R helper can load both stores and produce expected data frames/matrices.

Do not add downstream tests for `--save-intermediates` until those files become
part of a deliberate public analysis contract.

## Documentation Updates

Update:

- CLI help in `ref_gc_bias/config.rs`
- CLI help in `gc_bias/config.rs`
- shared `--gc-file` help in `cli_common.rs`
- generated CLI docs
- `website/docs/get-started/common-files.md`
- `.AI/docs/specs/gc_bias_flow_spec.md`
- command pipeline docs for `ref-gc-bias` and `gc-bias`
- py-cfdnalab README after helper support exists

## Implementation Order

1. Implement Zarr writing and reading for `ref-gc-bias`.
2. Change `ref-gc-bias` final output to `.ref_gc_package.zarr`.
3. Update `gc-bias --ref-gc-file` loading and validation to require the Zarr
   reference package.
4. Implement Zarr writing and reading for `GCCorrectionPackage`.
5. Change `gc-bias` final output to `.gc_bias_correction.zarr`.
6. Update downstream `--gc-file` loading and validation to require the Zarr
   correction package.
7. Update docs and generated CLI help.
8. Add Python helper support for both package types.
9. Add downstream compatibility checks for the two package stores.
10. Distill the final stable schema into `.AI/docs/specs/gc_bias_flow_spec.md`.

## Open Decisions

- Whether to keep the shared `GC_CORRECTION_SCHEMA_VERSION` for both Zarr
  package schemas or split into `REFERENCE_GC_PACKAGE_SCHEMA_VERSION` and
  `GC_CORRECTION_PACKAGE_SCHEMA_VERSION`. Prefer split constants because the
  two packages can evolve independently.
- Whether to write a one-time converter from `.npz` to `.zarr`. Do not add it
  unless backwards compatibility is explicitly needed.
- Whether root attrs should include the full command version once the project
  has a stable release-version source.
