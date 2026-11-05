# py-cfdnalab Loader Review

Date: 2026-11-05

Scope: `py-cfdnalab/src/cfdnalab/*`, the unit tests in `py-cfdnalab/tests/`,
`downstream_tests/python/`, `py-cfdnalab/README.md`, and the generated website
loader pages under `website/docs/generated/loaders/`. Claims about output
semantics were checked against the Rust writers.

Method: everything was read. Claims that depend on library behavior or on
performance at realistic sizes were checked with small probe scripts in the
package venv (pandas 3.0.3, zarr 3.1.6, numpy 2.4.5). The test suite was not
run. Anything marked *inference* was not measured.

Finding numbers match the first chat version of this review, so earlier
discussion can be traced. Items you rejected in chat were removed or rewritten.

## Checked And Not Findings

- Every hand-derived expectation in `tests/test_reference_correction.py` and
  `tests/test_reference_correction_core.py` was re-derived by hand and is
  correct.
- The correction formula `count / (reference_frequency * positive_support)`
  matches the Rust loader. Split mode multiplies the two side denominators.
  Outside and inside modes aggregate counts by side first.
- Canonical k-mer labels match Rust `make_canonical` (`src/shared/base.rs:130`)
  for even and odd k. Odd k is covered by the k=3 canonical downstream fixture.
- Motifs-file reference frequencies are normalized over the motifs-file
  targets, which the correction code expects.
- Group index order is deterministic from the BED file.
- Long-format lengths frames are fast enough. 290k windows x 50 bins:
  `read_lengths` 1.1 s, wide frame 0.3 s, long frame 6.0 s.
- `use_global_bias` naming is fine. A global reference composition is the
  global bias.

## Release Triage

Fix before validation:

- P-01: `read_lengths` fails when every group name is numeric.
- P-02: The README correction example fails on default `cfdna ends` output.
- P-03: Sparse motif axes differ between samples. Needs a design decision.
- P-04 and P-05: Midpoint values and end-motif labels are not explained.

Fix before people use it on genome-wide data:

- P-06 to P-08: Midpoint profile reads, sparse data frames, and reference
  correction get slow at realistic sizes.

## Wrong Results Or Broken Use

### P-01 - High - `read_lengths` fails when every group name is numeric

Evidence:

- `py-cfdnalab/src/cfdnalab/lengths.py:786` calls `pd.read_csv` without
  column types, so pandas guesses them.
- `py-cfdnalab/src/cfdnalab/lengths.py:844-850` converts only `chrom` back to
  text.
- `py-cfdnalab/src/cfdnalab/lengths.py:1033-1040` then rejects anything that
  is not a `str`.

Impact:

A grouped BED whose name column holds IDs like `1`, `2` gives integer group
names, and loading fails with
`ValueError: group_name must contain non-empty character strings` (probe
confirmed). Names like `001` only survive when the file also has a non-numeric
name.

Fix:

Pass `dtype={"chrom": str, "group_name": str}` to `read_csv` and remove
`_normalize_length_count_table`. Add a test with all-numeric group names. The
current test only covers numeric chromosomes.

### P-02 - High - README correction example fails on default `cfdna ends` output

Evidence:

- `cfdna ends` writes sparse output unless `--all-motifs` is set
  (`src/commands/ends/ends.rs:651`).
- `py-cfdnalab/README.md:286-289` calls `corrected_counts_array(...)` without
  `allow_densify=True`.
- `py-cfdnalab/src/cfdnalab/reference_correction.py:427-432` refuses to turn
  sparse output dense without that flag.

The website guide works only because it runs `cfdna ends --all-motifs`.

Fix:

Add `allow_densify=True` to the example, or show
`sparse_corrected_counts_matrix` instead.

### P-03 - High - Sparse motif axes differ between samples

Evidence:

- In sparse end-motif and ref-kmer output, the motif axis holds only the
  motifs observed in that file. It is sorted for full k-mer runs and follows
  the motifs-file order for `--motifs-file` runs
  (`downstream_tests/R/test_cfdnalab_r_package.R:230` expects
  `GT_AC`, `AC_GT`).
- Requesting an unobserved motif with `motifs=[...]` raises `KeyError`.
- `motif_index` and array column positions are local to each file.
- Only `website/docs/guides/fragment_ends_guide.md:37-39` warns about this.
  The README and docstrings do not.

Impact:

Comparing samples is the main use of end motifs, and nothing in the package
helps line samples up:

- Two samples can return arrays with the same number of columns but different
  motif sets. Stacking them silently misaligns columns.
- A user who already has a motif set of interest (a published panel, or the
  motif list from a trained model) cannot ask a sample for it if some motifs
  were not observed in that sample.
- Long data frames are safe when joined on `motif`, but not on
  `motif_index`. Nothing says which to use.

Design notes:

Expanding to every possible motif is not an option. 25-mers have about
`10^15` possible labels. Even short motifs on both sides multiply fast. Any
fix must keep the motif set explicit and bounded. Two separate pieces seem
needed:

1. A helper that finds a motif set across samples before anything is
   densified. It only needs each file's motif labels, not the counts, so it is
   cheap to run over a cohort. Useful rules:
   - Union: every motif observed in any sample. Missing entries become zero.
     For long motifs the union can be millions of labels, so the result must
     stay sparse (shared columns, sparse matrix). Densifying stays opt-in.
   - Intersection: motifs observed in every sample. Small, but it changes when
     a sample is added, so it is a poor fit for train/apply workflows.
   - Observed in at least `n` samples (or a fraction). Cuts the long tail of
     singletons from long motifs. A count-based threshold would need to read
     counts, which costs more.

2. Loader methods that accept an explicit motif list and fill motifs absent
   from this file. This covers both a set from the helper above and a set the
   user already has. Without it, the helper is not enough, because a model
   trained on one cohort must be applied to new samples with the same columns.

Filling a missing motif with zero is only correct when the run could have
counted that motif. The loader can check part of this and not the rest:

- Label shape: requested labels must match this file's outside and inside
  lengths and use valid bases. The loader can read the shape from the stored
  labels.
- Motifs-file runs: a motif outside the motifs file was never counted, so it
  is unknown, not zero. The Zarr root does not record the target list. Only
  `end_settings.json` records the motifs-file path, not its contents
  (`src/commands/ends/write.rs:53-62`). Either record the target list (or a
  hash of it) in the store, or make the user confirm the list.
- Motif-group axes: group names come from the user's group file. Same problem
  as motifs-file runs.
- Hidden `--collapse-complement` runs: a non-collapsed label would be filled
  with zero when it was really merged into its complement. See P-27.

An option name like `missing_motifs="zero"` (default `"error"`) keeps today's
behavior unless the user asks for filling.

Also document now, independent of the above, that `motif_index` is local to a
file and that samples should be joined on `motif`.

### P-04 - High - Midpoint values and positions are not explained

Nothing in the README, the docstrings, or the website pages generated from them
says what the `count` column holds or where positions start. The values are:

- Weighted (GC and scaling) `f32` values (`src/commands/midpoints/zarr.rs:128`).
- Smoothed by default with an order-3 Savitzky-Golay filter over 165 bp
  (`src/commands/midpoints/config.rs:176`), so they can be negative.
- Per-bp averages when `--bin-size > 1`
  (`src/commands/midpoints/postprocess.rs:197`), not sums.
- Summed over the intervals in a group. Dividing by `eligible_intervals`
  gives the mean per site.

Positions are 0-based offsets from the interval start, mirrored for `-` strand
intervals (`src/commands/midpoints/strand.rs:16`). No coordinate relative to
the interval center is offered, which is what most users plot.

The unit-test fixture adds to the confusion. It uses center-relative positions
`[-4, -2, 0, 2]` (`py-cfdnalab/tests/test_midpoints.py:20-21`), which the
writer never produces (`src/commands/midpoints/zarr.rs:262`). The downstream
test confirms real output starts at 0 (`[0, 2, 4, 6, 8]`).

Fix:

- Document the points above in `read_midpoints`, `data_frame`, and the README.
- Consider a center-relative position column in `positions()` and
  `data_frame()`.
- Consider a `settings()` method that returns `midpoint_settings.json`. It
  already records smoothing, GC correction, and scaling. Note that
  `r-cfdnalab/R/lengths.R:26-28` records a decision to keep settings JSON out
  of the lengths loader "so command provenance stays separate from the tabular
  data contract". If that decision covers every loader, drop this idea and let
  the docstrings carry the explanation instead.

### P-05 - High - End-motif labels and storage are not explained

The Python docs do not tell users:

- What `AC_GT` means: 2 outside bases from the reference, then 2 inside bases
  (from the read by default), read 5' to 3' inward from each fragment end.
  Right-end motifs are reverse-complemented.
- That sparse output holds observed motifs only and dense output means the run
  used `--all-motifs`. The README says to "start by checking storage mode"
  without saying why the two differ.
- That counts are weighted floats, not integers.

## Performance At Realistic Sizes

### P-06 - High - Getting midpoint profiles for many groups is slow

Evidence:

`py-cfdnalab/src/cfdnalab/midpoints.py:272-276` fetches one profile at a time:
one request for each selected group and length bin. `data_frame()` uses the
same path.

What happens:

The Zarr store keeps counts in large compressed blocks. The Rust writer
targets about 4 million values per block, covering many groups. Each
request unpacks the whole block, keeps the one profile it asked for, and
throws the rest away. Asking for 300 groups x 20 length bins unpacks the same
few blocks 6000 times.

Measured on a 300-group x 20-length-bin x 2001-position store with the writer's
block layout:

- 10 groups: 3.5 s with the current code, 0.02 s with one request for all of
  them.
- All 300 groups: about 100 s with the current code (scaled from the 10-group
  run).

For a user, this is the "plot every TSS group" notebook cell taking minutes,
and the conclusion that the package is slow.

Fix:

Ask Zarr once for all selected groups and length bins together, so each block
is unpacked once. In zarr-python this is
`counts.oindex[group_indices, length_bin_indices, :]`. The public API does not
change.

### P-07 - Medium - Sparse data frames scan the whole file in Python, even for one window

Evidence:

`py-cfdnalab/src/cfdnalab/ends.py:762-781` loops in Python over every stored
entry to find the requested rows and motifs.
`py-cfdnalab/src/cfdnalab/ref_kmers.py:754-773` is a copy of it.

Measured: `data_frame(window_idxs=0)` took 0.56 s on a store with 5.1M stored
entries. The cost grows with the size of the file, not the selection, so a
loop over windows becomes quadratic.

Fix:

Replace the loop with array lookups: a table from stored row to selected
position, then a mask and a sort.

### P-08 - High - Reference correction is slow and memory-hungry on genome-wide windows

Evidence:

`py-cfdnalab/src/cfdnalab/reference_correction.py:284-396` always builds long
pandas tables and joins them on text keys. The array and matrix outputs go
through the same path (`reference_correction.py:399-550`).

Measured with 256 motifs and full support:

- 3k windows: 2.2 s.
- 20k windows: 20 s and about 1.5 GiB extra memory.

*Inference*: genome-wide 10 kb windows (about 290k rows, the window size the
README recommends) would take about 5 minutes and likely more than 16 GiB.

Fix:

Do it with matrices, like the Rust loader:

- Line the reference frequencies up with the sample's motif columns once.
- Count supported motifs per row (`frequency > 0`).
- Divide the counts (or the stored sparse values) by
  `frequency * supported_motif_count`.
- Get side frequencies for split, outside, and inside modes by summing
  reference columns that share a side label.

The data frame output can be built from the result.

### P-09 - Low - Smaller costs

- Dense `sparse_counts_matrix(window_idxs=...)` reads the whole dense matrix
  before selecting rows (`py-cfdnalab/src/cfdnalab/ends.py:343-345`).
- Loading sparse ref-kmers checks entry order in a Python loop
  (`py-cfdnalab/src/cfdnalab/ref_kmers.py:1547-1554`). 0.9 s per 5M entries.
- Neither sparse data frame helper keeps the sparse matrix it builds, so
  repeated calls rebuild it.

## API Consistency

### P-10 - Medium - Index column names are inconsistent

Data frames use `motif_index`, `length_bin`, and `position`. Selectors use
`motif_idx()`, `motif_idxs`, and `length_bin_idxs`. Meanwhile the
`group_idx` and `window_idx` columns do use `_idx`. A user writing
`meta["motif_idx"]` gets a `KeyError`.

The R package already uses `motif_idx`, `length_bin_idx`, and `position_idx`
columns, so the two packages also disagree with each other. A user moving
between them meets different column names for the same thing.

Suggestion: use `*_idx` for every index column.

### P-11 - Low - `length_bins()` returns different columns for lengths and midpoints

What each returns:

- `LengthCounts.length_bins()` (`py-cfdnalab/src/cfdnalab/lengths.py:112-137`):
  `length_bin`, `length_start_bp`, `length_end_bp`, `length_midpoint_bp`,
  `length_width_bp`.
- `MidpointProfiles.length_bins()`
  (`py-cfdnalab/src/cfdnalab/midpoints.py:279-297`): `length_bin`,
  `length_start_bp`, `length_end_bp`.

Why it may be on purpose:

In lengths, the bins are the x-axis of a length distribution. The midpoint is
where a bar or point is drawn, and the width is needed for `value="density"`
when bins have unequal widths. In midpoints, length bins are usually a few
broad classes (for example short and long fragments) used to pick which
profile to look at. There, midpoint and width add little.

Why it may still be worth adding:

A midpoint run with many narrow length bins gives a fragment length x position
heatmap (a V-plot). For that, users need a y-coordinate per length bin, which
is `length_midpoint_bp`. They can compute it themselves, but lengths already
gives it to them.

Question for you: keep the difference, add the two columns to midpoints, or
mention the difference in the midpoints docstring.

### P-12 - Medium - Corrected arrays filtered by `max_blacklisted_fraction` do not say which rows they hold

Evidence:

- The windowed and grouped `corrected_counts_array` and
  `sparse_corrected_counts_matrix` accept `max_blacklisted_fraction`
  (`py-cfdnalab/src/cfdnalab/ends.py:1569`, `1676`, and the grouped versions).
- Rows above the cutoff are dropped
  (`py-cfdnalab/src/cfdnalab/reference_correction.py:441-447`, `500-506`).
- The docstrings do not mention the option (see P-16).
- The returned array has no row labels. The raw `dense_counts_array` does not
  take this option, and `window_metadata()` / `group_metadata()` take no
  arguments.

Impact:

The filter is useful and should stay. The problem is that after filtering, a
user cannot easily tell which windows the rows are, and cannot get a raw array
with the same rows. Plotting corrected values along the genome then means
re-implementing the filter (`blacklisted_fraction <= cutoff`, applied after
the window selection).

Options that keep the feature:

- Document the option and the row-dropping behavior on the corrected methods.
- Let `window_metadata()` and `group_metadata()` take the same row selector and
  `max_blacklisted_fraction`, so the same arguments give matching metadata.
  This is small and fits the existing API.
- Let the raw array methods take `max_blacklisted_fraction` too, so raw and
  corrected arrays can be filtered the same way.

### P-13 - Discussion - Raw end-motif data frames have no frequency column

You mentioned there was a reason. I found no recorded reason in the specs, the
CHANGELOG, or the git history of `ends.py`. Facts that matter for the
discussion:

- Ref-kmer frames have `frequency`. Corrected frames have
  `corrected_frequency`.
- `corrected_frequency` is normalized over the full stored motif axis before
  the motif selector is applied
  (`py-cfdnalab/src/cfdnalab/reference_correction.py:375-382`). So selecting
  motifs does not change the reported frequencies. The same rule could serve
  raw counts.
- For sparse output, the full stored axis is every observed motif, so the row
  total is still the total over all counted ends.
- *Inference* on possible reasons: for motifs-file runs, a frequency is
  relative to the targeted motifs, not to all fragment ends. For motif-group
  axes it is relative to the groups. A user could misread either as a share of
  all ends.

### P-14 - Low - Unknown labels raise different error types

An unknown motif raises `KeyError`, but an unknown side label in side mode
raises `ValueError` (`py-cfdnalab/src/cfdnalab/reference_correction.py:774`).

### P-15 - Low - `row_mode()` values follow CLI flag names

`EndMotifCounts.row_mode()` and `RefKmerFrequencies.row_mode()` return
`"size"`, `"bed"`, or `"grouped_bed"`. These mirror `--by-size`, `--by-bed`,
and `--by-grouped-bed` (`src/commands/cli_common.rs:408-437`,
`src/commands/ends/zarr.rs:287-292`). The classes are named `Windowed...` and
`Grouped...`. Lengths has no `row_mode()` and its TSV does not record size
versus BED windows.

This is probably fine, since `"size"` versus `"bed"` carries information that
`"windowed"` would lose. At most, say in the `row_mode()` docstring which flag
each value comes from and which class it gives.

## Documentation

### P-16 - Medium - Corrected-count methods and all ref-kmer methods lack Parameters and Returns

This is not docstring inheritance. Checked:

- Each of the six `corrected_counts_array` and `sparse_corrected_counts_matrix`
  methods defines its own docstring
  (`py-cfdnalab/src/cfdnalab/ends.py:1126`, `1233`, `1561`, `1669`, `2034`,
  `2143`). None has a Parameters or Returns section. `allow_densify`,
  `max_blacklisted_fraction`, `use_global_bias`, `unsupported_motifs`, and
  `two_sided_correction` are undocumented there.
- All 32 public ref-kmer methods define their own docstrings and none has
  Parameters or Returns. That includes the ones that take arguments, such as
  `data_frame`, `motif_idx`, `group_idx`, and the array and matrix methods.
  `read_ref_kmers`
  (`py-cfdnalab/src/cfdnalab/ref_kmers.py:1287`) has a 3-line docstring.
  `read_end_motifs`, `read_lengths`, and `read_midpoints` have both sections.
- The website generator parses each class body statically with
  `ast.get_docstring` (`website/scripts/generate_loader_docs.py:789-829`), so
  nothing would be inherited there either. The website pages have the same
  gaps.

Suggestion: `validate_python_docstrings`
(`website/scripts/generate_loader_docs.py:1301-1327`) only checks that a
docstring exists. It could also require a Parameters section when the method
takes arguments, so this cannot come back.

### P-17 - Low - The reference-correction explanation is copied into 9 docstrings

The about 60-line "Reference correction" and "Two-sided correction modes"
sections appear in 9 docstrings in `ends.py`. The copies are identical today.

Ways to keep them from drifting:

- Keep the copies and add a test that compares those sections across the 9
  docstrings. No change for readers. Cheapest.
- Write the full explanation once (for example in the website ends guide) and
  replace each copy with a short summary and a pointer. Shorter docstrings, but
  `help(...)` no longer shows the full explanation.
- Build the docstrings at import time from one shared string. `help(...)`
  would work, but the website generator reads source statically and would show
  the unexpanded text. Not recommended without generator changes.

I would start with the test.

### P-18 - Notes - README

Notes only:

- "was used" should be "were used". "why we recommend" should be "which is why
  we recommend".
- The eligible-intervals example (`py-cfdnalab/README.md:78-88`) loops group by
  group with `iterrows`. One call with `group_idxs=` does the same.
- End-to-end recipes would help: a TSS profile normalized to its flanks, a
  short/long ratio per window, a samples x motifs frequency matrix.

### P-19 - Low - Length counts are called "raw" but can be weighted

`counts_array` says "raw length counts". They are weighted when GC correction
or scaling is on. `length_settings.json` and `end_settings.json` record this,
and the loader does not surface either file (see the `settings()` idea in
P-04).

## Tests

- P-20: The midpoint fixture uses positions the writer never produces (P-04).
- P-21: No test covers all-numeric group names (P-01).
- P-22: Corrected values for motif-group outputs are not tested anywhere: not
  in the Python unit tests, the downstream tests, or the Rust loader tests.
  Only the rejection of two-sided modes is tested.
- P-23: Also untested:
  - Corrected arrays with `max_blacklisted_fraction` (rows dropped).
  - `keep_na` NaN values stored in the sparse corrected matrix.
  - `corrected_counts_array` on sparse input with `allow_densify=True`.
- P-24: Unit fixtures mostly use end-motif schema version 1 without
  `motif_axis_kind`. The writer produces version 2. The downstream tests cover
  version 2.
- `tests/test_package.py` only checks that public symbols exist.
- The midpoint group name `"gamma-unicode-aa"` suggests a Unicode test, but the
  name is ASCII.

## Design Questions

### P-25 - Low - A sparse end-motif store with zero counts cannot be opened

`py-cfdnalab/src/cfdnalab/ends.py:2308-2327` rejects it, so even the window
metadata is unreadable. A dense store with zero counts opens fine.

The R loader has to reject these stores, because the CRAN `zarr` reader cannot
open zero-length arrays (`r-cfdnalab/R/end_motifs.R:226-235`). If the Python
rejection exists to match R, a comment saying so would help. If not, Python
could open the store and return empty results.

### P-26 - Low - Blacklist cutoffs can keep different windows in lengths and ends

Lengths TSVs round `blacklisted_fraction` to 3 decimals
(`src/commands/lengths/writer.rs:24`). Zarr outputs keep full precision. A
window with `0.1004` blacklisted is `0.100` in lengths and passes a `0.1`
cutoff there, but fails it in ends.

### P-27 - Low - `--collapse-complement` is not recorded in the Zarr store

The hidden flag is only written to `end_settings.json`
(`src/commands/ends/write.rs:127`). Correction cannot detect it and would apply
single k-mer frequencies to collapsed labels.

### P-28 - Low - Some errors are hidden or misreported

- `_has_array` (`py-cfdnalab/src/cfdnalab/ends.py:2446-2454`) swallows every
  exception, so a damaged store is reported as a missing array.
- The `zstd` command-line fallback in `_read_length_counts_table`
  (`py-cfdnalab/src/cfdnalab/lengths.py:777-822`) can never run, because
  `zstandard` is a required dependency. Its broad `except` also hides the real
  parse error for a malformed `.zst`.

### P-29 - Low - Loading code is duplicated and validation has drifted apart

`ends.py`, `ref_kmers.py`, and `midpoints.py` each carry their own loading and
selector code. `ref_kmers` checks dimensions, dtypes, intervals, fraction
ranges, and sparse entry order. `ends` checks much less. A shared module would
fix this, but it is a large refactor and needs your go-ahead.
