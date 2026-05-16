# Python Package Plan

This plan covers a small Python companion package for reading and analyzing
cfDNAlab output files. The first useful target is midpoint Zarr output.

## Purpose

The package should make public cfDNAlab outputs convenient to use from Python
after the Rust CLI has already produced them. It is not a Python implementation
of cfDNAlab and should not pretend to install or wrap the Rust command-line
tool.

The initial package should focus on:

- opening `<prefix>.midpoint_profiles.zarr`
- validating the midpoint Zarr schema and version
- exposing compact group, length-bin, and position metadata
- extracting NumPy arrays for common slices
- building pandas data frames for common plotting and analysis shapes

Do not add plotting in the first package version. The first useful boundary is
clean arrays and data frames.

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
py_cfdnalab/
  pyproject.toml
  README.md
  src/
    cfdnalab/
      __init__.py
      midpoints.py
  tests/
    test_midpoints.py
```

Use `src/` layout so local imports behave like an installed package and do not
accidentally import from the repository root.

The Rust crate should not depend on the Python package. The Python package can
live in the same repository, but it should have its own metadata, tests, and
development commands.

## Dependencies

Initial runtime dependencies:

```toml
dependencies = [
  "numpy",
  "pandas",
  "zarr>=3",
]
```

Use `requires-python = ">=3.11"` because Zarr Python 3 requires Python 3.11
or later.

Initial test dependency:

```toml
[project.optional-dependencies]
test = [
  "pytest",
]
```

Do not add xarray, dask, matplotlib, seaborn, or plotnine as runtime
dependencies for the first version. Those can remain user-side plotting choices
or downstream compatibility checks.

## Public API

The package-level API should be:

```python
import cfdnalab as cfl

midpoints = cfl.load_midpoints("sample.midpoint_profiles.zarr")
```

`cfdnalab.__init__` should re-export only stable public entry points:

```python
from .midpoints import MidpointProfiles, load_midpoints
```

Initial midpoint helper methods:

```python
midpoints.groups()
midpoints.group_names()
midpoints.eligible_intervals()
midpoints.length_bins()
midpoints.positions()

midpoints.group_idx(group_name="CTCF")
midpoints.length_bin_idx(length=167)

midpoints.data_frame_for_profile(group_idx=0, length_bin=0)
midpoints.data_frame_from_group(group_name="CTCF")
midpoints.data_frame_from_group_idx(group_idx=0)
midpoints.data_frame_from_length(length=167)
midpoints.data_frame_from_length_bin(length_bin=0)

midpoints.array_for_profile(group_idx=0, length_bin=0)
midpoints.array()
midpoints.array_from_group(group_name="CTCF")
midpoints.array_from_group_idx(group_idx=0)
midpoints.array_from_length(length=167)
midpoints.array_from_length_bin(length_bin=0)
```

`array()` must be documented as loading the full 3D tensor into RAM. Examples
should prefer `array_for_profile`, `array_from_group`, or
`array_from_length_bin` when possible.

## Tests

Use two test layers.

Python package unit tests:

- live under `py_cfdnalab/tests`
- create a tiny valid Zarr store directly, or use a small checked-in fixture if
  direct fixture creation is less brittle
- test schema validation errors, metadata helpers, array slicing, dataframe
  builders, group-name lookup, and length-bin lookup
- run with:

```bash
cd py_cfdnalab
python -m pip install -e ".[test]"
pytest
```

Downstream integration tests:

- keep using cfDNAlab-generated Zarr output
- verify that the Python package can read real Rust CLI output
- should run in the existing downstream workflow or a Python-package-specific
  job that first generates the Rust fixture

The package unit tests should not require BAM generation or the Rust CLI. The
integration tests are the layer that catches Rust output-schema drift.

## README Requirements

The full package README should include:

- what this package is
- what it is not
- how to install the Rust CLI separately
- how to install the Python package
- how to load midpoint Zarr output
- how to inspect group, length-bin, and position metadata
- how to get one profile as a data frame
- how to filter groups by `eligible_intervals`
- how to get NumPy arrays
- a warning that `array()` loads the full 3D count tensor into RAM

The first example should show the split between producing output and reading it:

```bash
cfdna midpoints ...
pip install cfdnalab
```

```python
import cfdnalab as cfl

midpoints = cfl.load_midpoints("sample.midpoint_profiles.zarr")
profile = midpoints.data_frame_for_profile(group_idx=0, length_bin=0)
```

## Release Notes

When this package is published, add a note to the main cfDNAlab documentation:

- PyPI `cfdnalab` is a Python reader/helper package
- the Rust CLI is still installed/distributed separately
- the Python package is intended for downstream analysis of output files

## Later Extensions

Possible later additions:

- R package with equivalent midpoint helpers
- ends Zarr reader helpers after the ends schema is implemented
- optional plotting helpers if common plot shapes stabilize
- xarray export method if it proves useful in real workflows

Do not add these before the midpoint reader API has been tested on real outputs.
