# R Package Development

This document is for developers working on the `r-cfdnalab` package. The public
README should stay focused on installation and ordinary use.

## Documentation

Documentation is written with roxygen comments in `R/*.R`. Regenerate
`NAMESPACE` and `man/*.Rd` with:

```r
roxygen2::roxygenise("r-cfdnalab")
```

Do not edit generated `.Rd` files directly. Update the roxygen comments instead
and regenerate the generated files.

## Indexing

The Zarr schemas store axis coordinates as zero-based integer arrays because
the files are shared with Python and other non-R readers. The R package must
not expose those zero-based coordinates as ordinary user selectors.

Public R selectors and public metadata data frames are one-based:

- Use public names like `group_idx`, `length_bin_idx`, `position_idx`,
  `window_idx`, `chromosome_idx`, `motif_idx`, and `row_idx`.
- `group_idx(x, "alpha")`, `length_bin_idx(x, 167)`, and `motif_idx(x, "_A")`
  return one-based R indices.
- Public data frames from helpers such as `groups()`, `length_bins()`,
  `positions()`, `windows()`, `motifs()`, and sparse/dense data-frame helpers
  must contain one-based `*_idx` columns.
- Public data frames must not contain internal `*_idx0` or `*_index0` columns.

Internal zero-based values must be explicit in the name:

- Store raw schema axes as fields such as `group_idx0`, `length_bin_idx0`,
  `position_idx0`, `row_idx0`, and `motif_idx0`.
- Use helper-local variables with the same convention, for example
  `group_idx0` or `motif_idx0`.
- Keep schema array names unchanged. For example, the end-motif Zarr array is
  still named `motif_index`; only the R public method is `motif_idx()`.

Use the shared conversion helpers instead of open-coded conversions:

```r
cf_validate_r_index(index, size, name)
cf_r_index_to_index0(index)
cf_index0_to_r_index(index0)
cf_validate_index0(index0, size, name)
```

Tests should include asymmetric cases where public index `1` and internal
zero-based index `1` select different values. That catches accidental exposure
of schema indices as R-facing indices.
