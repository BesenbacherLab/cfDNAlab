# Python Package Plan

This plan covers the Python companion package for reading and analyzing
cfDNAlab output files. The current public loader API is specified in
[`../specs/package_loader_api.md`](../specs/package_loader_api.md).

## Purpose

The package should make public cfDNAlab outputs convenient to use from Python
after the Rust CLI has already produced them. It is not a Python implementation
of cfDNAlab and should not pretend to install or wrap the Rust command-line
tool.

The package currently targets:

- `<prefix>.midpoint_profiles.zarr`
- `<prefix>.end_motifs.zarr`
- `<prefix>.length_counts.tsv.zst`

The useful boundary is clean metadata, NumPy arrays, SciPy sparse matrices, and
pandas data frames. Plotting can stay outside the first package scope until
common display shapes stabilize.

## Naming

Use the PyPI distribution name `cfdnalab` if available, to reserve the project
name and avoid third-party squatting. Use the Python import name `cfdnalab`.

The package description must be explicit:

```text
Python helpers for loading cfDNAlab output files
```

Avoid descriptions like "cfDNAlab for Python", because that implies the Rust CLI
is installed through pip.

The README must state near the top:

```text
This package does not install the cfDNAlab command-line tool. The CLI is
distributed separately as the Rust `cfdna` binary. Use this Python package after
running cfDNAlab to load and analyze output files.
```

Do not add console scripts named `cfdna` or `cfdnalab` in the Python package.
That would increase confusion with the Rust CLI.

## Conda And Bioconda Naming

Keep conda package names explicit about which runtime they provide.

Preferred future split:

```text
cfdnalab      Rust CLI package; installs the `cfdna` executable
py-cfdnalab   Python helper package; imports as `cfdnalab`
r-cfdnalab    R helper package
```

`cfdnalab` on Bioconda should mean the main Rust command-line tool, not the
Python helper. This matches user expectations for a project-level package name
and keeps CLI-only environments free of Python and R analysis dependencies.

The Python and R helper packages should depend on their own language ecosystem
dependencies only. They should read cfDNAlab output files; they should not wrap
or vendor the Rust CLI.

If users eventually want everything in one environment, add a separate
metapackage instead of making the core CLI package heavy:

```text
cfdnalab-suite
```

That metapackage could depend on `cfdnalab`, `py-cfdnalab`, and `r-cfdnalab`.
Do not make this the default package unless there is clear user demand.

## Repository Layout

Keep the Python package as an isolated subproject:

```text
py-cfdnalab/
  pyproject.toml
  README.md
  src/
    cfdnalab/
      __init__.py
      _helpers.py
      ends.py
      lengths.py
      midpoints.py
  tests/
    test_ends.py
    test_lengths.py
    test_midpoints.py
```

Use `src/` layout so local imports behave like an installed package and do not
accidentally import from the repository root.

The Rust crate should not depend on the Python package. The Python package can
live in the same repository, but it should have its own metadata, tests, and
development commands.

## Dependencies

Runtime dependencies:

```toml
dependencies = [
  "numpy",
  "pandas",
  "scipy",
  "zarr>=3",
  "zstandard",
]
```

Use `requires-python = ">=3.11"` because Zarr Python 3 requires Python 3.11 or
later.

Test dependency:

```toml
[project.optional-dependencies]
test = [
  "pytest",
]
```

Do not add xarray, dask, matplotlib, seaborn, or plotnine as runtime
dependencies. Those can remain user-side plotting choices or downstream
compatibility checks.

## Public API

`cfdnalab.__init__` should re-export stable public entry points:

```python
from .ends import (
    EndMotifCounts,
    GlobalEndMotifCounts,
    GroupedEndMotifCounts,
    WindowedEndMotifCounts,
    read_end_motifs,
)
from .lengths import (
    GlobalLengthCounts,
    GroupedLengthCounts,
    LengthCounts,
    WindowedLengthCounts,
    read_lengths,
)
from .midpoints import MidpointProfiles, read_midpoints
```

The package-level API should stay simple:

```python
import cfdnalab as cfl

midpoints = cfl.read_midpoints("sample.midpoint_profiles.zarr")
ends = cfl.read_end_motifs("sample.end_motifs.zarr")
lengths = cfl.read_lengths("sample.length_counts.tsv.zst")
```

The mode-specific data frame methods and selector names are specified in
[`../specs/package_loader_api.md`](../specs/package_loader_api.md). Keep Python
as object-method-oriented: `.data_frame(...)` is appropriate because it is not a
global function that risks masking another package.

Array and matrix helpers should make large or dense materialization visible:

- midpoint `.counts_array()` returns a 3D NumPy array with group, length-bin,
  and position axes
- end-motif `dense_counts_array()` may load or reconstruct dense count arrays
- end-motif `sparse_counts_matrix()` preserves sparse output as a SciPy sparse
  matrix
- length `.counts_array()` returns a 2D NumPy array with output-row and
  length-bin axes

## Tests

Use two test layers.

Python package unit tests:

- live under `py-cfdnalab/tests`
- create small valid package fixtures directly, or use small checked-in package
  fixtures if direct construction is less brittle
- test schema validation errors, metadata helpers, array/matrix extraction,
  data frame builders, selector validation, compressed length TSV reading, and
  duplicate-name handling
- run with:

```bash
cd py-cfdnalab
python -m pip install -e ".[test]"
pytest
```

Downstream integration tests:

- use cfDNAlab-generated midpoint, end-motif, and length fixtures
- verify that the Python package can read real Rust CLI output
- include fixture variants where no blacklist was used, so loader behavior does
  not accidentally require `blacklisted_fraction`
- should run in the existing downstream workflow or a package-specific job that
  first generates the Rust fixture

The package unit tests should not require BAM generation or the Rust CLI. The
integration tests are the layer that catches Rust output-schema drift.

## README Requirements

The package README should include:

- what this package is
- what it is not
- how to install the Rust CLI separately
- how to install the Python package
- how to load midpoint, end-motif, and length-count outputs
- how to inspect group, length-bin, window, position, and motif metadata
- how to get representative data frames without documenting every method
- how to get NumPy arrays and SciPy sparse matrices
- a warning that full-array helpers can load large data into RAM
- a note that stricter `max_blacklisted_fraction` cutoffs require
  `blacklisted_fraction` metadata in the loaded output

The first example should show the split between producing output and reading it:

```bash
cfdna midpoints ...
pip install cfdnalab
```

```python
import cfdnalab as cfl

midpoints = cfl.read_midpoints("sample.midpoint_profiles.zarr")
profile = midpoints.data_frame(group_idxs=0, length_bin_idxs=0)
```

Keep examples representative, not exhaustive. The README should show the main
metadata, array, sparse-matrix, and data frame paths without becoming a full API
reference.

## Release Notes

When this package is published, add a note to the main cfDNAlab documentation:

- PyPI `cfdnalab` is a Python reader/helper package
- the Rust CLI is still installed/distributed separately
- the Python package is intended for downstream analysis of output files

## Later Extensions

Possible later additions:

- optional plotting helpers if common plot shapes stabilize
- xarray export method for midpoint profiles if it proves useful in real
  workflows
- GC-output loaders after those public schemas settle
