# py-cfdnalab Loader Review

Date: 2026-11-05

Scope: `py-cfdnalab/src/cfdnalab/*`, the unit tests in `py-cfdnalab/tests/`,
`downstream_tests/python/`, `py-cfdnalab/README.md`, the website guides, and
the generated website loader pages under `website/docs/generated/loaders/`.
Claims about output semantics were checked against the Rust writers.

Method: everything was read. Claims that depend on library behavior or on
performance at realistic sizes were checked with probe scripts in the package
venv (pandas 3.0.3, zarr 3.1.6, numpy 2.4.5). Probe stores used the Rust
writer's chunk layout and fill values. Memory figures are process peak resident
memory from `/usr/bin/time -l`, reported next to the peak after loading only.
The test suite was not run. Anything marked *inference* was not measured.

Finding numbers are stable across revisions of this review. Rejected items were
removed rather than renumbered.

## Checked And Not Findings

- Every hand-derived expectation in `tests/test_reference_correction.py` and
  `tests/test_reference_correction_core.py` was re-derived and is correct.
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
- `corrected_counts_array(..., allow_densify=True)` on sparse input is covered
  by writer-backed downstream tests for split and outside correction, with and
  without a motifs file
  (`downstream_tests/python/test_end_motif_zarr_cfdnalab_package.py:244-278`,
  `:374-408`).
- `corrected_motifs_metadata()` (`ends.py:598`) already gives the column
  identity of corrected arrays.

## Release Triage

Fix before validation:

- P-01: Numeric-looking chromosome and group labels are changed or rejected.
- P-02: README and website corrected-array examples fail on default output.
- P-03: Samples cannot be aligned on a shared motif axis. Document now, design
  the API afterwards.
- P-04 and P-05: Midpoint values and end-motif labels are not explained.
- P-27: Complement-collapsed output gets a wrong reference correction without
  an error. The guard is cheap.

Fix before people use it on genome-wide data:

- P-06 to P-08: Midpoint profile reads, sparse data frames, and reference
  correction get slow at realistic sizes.

## Wrong Results Or Broken Use

### P-01 - High - Numeric-looking chromosome and group labels are changed or rejected

Evidence:

- `py-cfdnalab/src/cfdnalab/lengths.py:786` calls `pd.read_csv` without
  column types, so pandas guesses them.
- `lengths.py:844-850` converts `chrom` back to text after parsing. By then
  the original text is gone.
- `lengths.py:1033-1040` rejects a `group_name` column that is not text.

Probes:

- Windowed chromosomes `001` and `1` both load as `"1"`. No error.
- Group names `1` and `2` fail with
  `ValueError: group_name must contain non-empty character strings`.
- Group names `1.10` and `1.1` fail the same way. Casting after parsing would
  have turned them into `"1.1"` and `"1.1"`.
- Names like `001` only survive when the same column also has a non-numeric
  value.

Impact:

Chromosome identity can be lost silently, and distinct labels can collapse into
one. Grouped outputs with numeric IDs in the BED name column cannot be loaded.

Fix:

Pass `dtype={"chrom": str, "group_name": str}` to `read_csv` and remove
`_normalize_length_count_table`. A cast after parsing must not be used for
either column, because it hides the error and keeps the lost identity.

Tests: all-numeric group names, zero-padded chromosomes (`001` and `1` in the
same file), and labels like `1.10` next to `1.1`. Cover plain and `.zst`
input and assert both the exact labels and selection by name. R-03 needs the
same cases.

### P-02 - High - README and website corrected-array examples fail on default output

Evidence:

- `cfdna ends` writes sparse output unless `--all-motifs` is set
  (`src/commands/ends/ends.rs:651`).
- `py-cfdnalab/README.md:286-289` and the Python block in
  `website/docs/guides/fragment_ends_guide.md:265-268` call
  `corrected_counts_array(...)` without `allow_densify=True`.
- `py-cfdnalab/src/cfdnalab/reference_correction.py:427-432` refuses to turn
  sparse output dense without that flag.
- The `--all-motifs` in the website guide (`fragment_ends_guide.md:217-227`)
  is passed to `cfdna ref-kmers`. It does not change the storage of the sample
  end-motif output.

Fix:

Show `sparse_corrected_counts_matrix()` in the general examples. Where a dense
array is shown, pass `allow_densify=True` with a bounded selection and say how
large the result gets.

### P-03 - High - Samples cannot be aligned on a shared motif axis

Evidence:

- In sparse end-motif and ref-kmer output, the motif axis holds only the
  motifs observed in that file. It is sorted for full k-mer runs and follows
  the motifs-file order for `--motifs-file` runs
  (`downstream_tests/R/test_cfdnalab_r_package.R:230` expects
  `GT_AC`, `AC_GT`).
- Requesting an unobserved motif with `motifs=[...]` raises `KeyError`.
- `motif_index` and array column positions are local to each file.
- The single-file facts are documented: `README.md:119` says densifying only
  includes observed motifs, and the `has_motif()` docstring explains absent
  motifs. What is missing is how to align several samples. Only
  `website/docs/guides/fragment_ends_guide.md:37-39` warns that motif sets
  differ between samples.
- Storage mode does not tell users whether the axis is complete. End-motif
  output is dense exactly when `--all-motifs` was used (`ends.rs:651`).
  Ref-kmer output is also dense when the observed axis happens to be complete
  (`src/commands/ref_kmers/ref_kmers.rs:488-489`).

Impact:

Comparing samples is the main use of end motifs, and nothing in the package
helps line samples up:

- Two samples can return arrays with the same number of columns but different
  motif sets. Stacking them misaligns columns without an error.
- A user with a fixed motif set (a published panel, or the motif list of a
  trained model) cannot request it from a sample that did not observe every
  motif.
- Long data frames are safe when joined on `motif`, but not on `motif_index`.
  Nothing says which to use.

Design notes:

Expanding to every possible motif is not an option. 25-mers have about
`10^15` possible labels. Any fix must keep the motif set explicit and bounded.
In order of priority:

1. Document alignment now: `motif_index` is local to a file, samples are joined
   on `motif`, and arrays must be reindexed to a shared motif order before
   stacking. State the completeness rule in terms of stored labels and the
   run's counting targets, not storage mode.
2. Let loader methods accept an explicit, bounded motif list and fill motifs
   absent from the file, behind an option such as `missing_motifs="zero"`
   (default `"error"`). This covers both a motif set the user already has and
   the train/apply workflow, where new samples need the training columns.
3. Later, a cohort helper that proposes a shared motif set (union,
   intersection, or observed in at least `n` samples). The union of long
   motifs can reach millions of labels, so results stay sparse. Intersection
   changes when a sample is added, which makes it a poor fit for train/apply.

Filling with zero is only correct when the run could have counted the motif.
What the loader can and cannot check:

- Label shape: requested labels must match the outside and inside widths and
  use valid bases. Widths can be read from stored labels, but an empty axis
  supplies none. Automatic validation needs the widths written to the store.
- Motifs-file runs: a motif outside the motifs file was never counted, so it
  is unknown, not zero. The Zarr root does not record the targets, and
  `end_settings.json` only records the motifs-file path
  (`src/commands/ends/write.rs:53-62`). The counting-target definition must be
  stored to validate membership. A hash can only confirm that two supplied
  definitions agree. It cannot recover an omitted target or decide membership.
- Motif-group axes: validation needs the motif-to-group mapping, not just
  matching group names.
- Complement-collapsed runs: see P-27.

Cost and comparability limits:

- A label-only cohort helper needs a metadata-only read path.
  `read_end_motifs()` and `read_ref_kmers()` read every sparse coordinate and
  value array when opening a store (`ends.py:162-168`, and the eager sorted
  check in `ref_kmers.py:1547-1554`).
- Ref-kmer frequencies are normalized over the run's counted targets. Zero
  filling cannot make runs with different targets comparable.

### P-04 - High - Midpoint values and positions are not explained

Nothing in the README, the docstrings, or the website pages generated from them
says what the `count` column holds or where positions start.

Values:

- Weighted `f32` values (`src/commands/midpoints/zarr.rs:128`). Optional GC and
  scaling weights are applied when fragments are counted.
- Smoothed by default with an order-3 Savitzky-Golay filter over 165 bp
  (`src/commands/midpoints/config.rs:176`), so values can be negative.
- Averaged per bp when `--bin-size > 1`, after smoothing
  (`src/commands/midpoints/postprocess.rs:193-215`). The last bin can be
  shorter and is divided by its real width.
- Summed over the retained intervals of each group. Dividing by
  `eligible_intervals` gives a mean per retained interval. Intervals can
  overlap or repeat, so this is not a mean over independent sites. Groups with
  no retained intervals stay in the output with `eligible_intervals = 0`
  (`src/commands/midpoints/group_index.rs:11-25`, `:79-83`), so divide only
  when it is positive.

Positions:

- 0-based offsets from the interval start.
- Mirrored for `-` strand intervals, where position 0 is the rightmost genomic
  base (`src/commands/midpoints/strand.rs:14-16`).
- No coordinate relative to the interval center is offered, which is what most
  users plot.

The unit-test fixture uses center-relative positions `[-4, -2, 0, 2]`
(`py-cfdnalab/tests/test_midpoints.py:20-21`), which the writer never produces
(`src/commands/midpoints/zarr.rs:262`). The downstream test confirms real
output starts at 0 (`[0, 2, 4, 6, 8]`). The fixture loads because the midpoint
loader does not validate bin intervals (P-29).

Fix:

- Document the points above in `read_midpoints`, `data_frame`, and the README.
  This does not require any API change.
- Optionally add a center-relative position column. It needs an explicit
  convention for even interval widths, for the point within a bin that the
  coordinate refers to, and for a shortened last bin.
- Optionally add a `settings()` method for `midpoint_settings.json`.
  `r-cfdnalab/R/lengths.R:26-28` records a decision to keep settings JSON out
  of the lengths loader "so command provenance stays separate from the tabular
  data contract". If that decision covers every loader, the docstrings must
  carry the explanation instead.

### P-05 - High - End-motif labels and storage are not explained

The Python docs do not explain:

- What `AC_GT` means: 2 outside bases from the reference, then 2 inside bases
  (from the read by default), read 5' to 3' inward from each fragment end.
  Right-end motifs are reverse-complemented.
- Why storage differs. End-motif output is dense exactly when `--all-motifs`
  was used. The README says to "start by checking storage mode" without saying
  why. For ref-kmers, dense storage does not imply `--all-motifs` (P-03).
- That counts are weighted floats, not integers.

## Performance At Realistic Sizes

### P-06 - High - Getting midpoint profiles for many groups is slow

Evidence:

`py-cfdnalab/src/cfdnalab/midpoints.py:272-276` fetches one profile at a time:
one request for each selected group and length bin. `data_frame()` uses the
same path.

What happens:

The Zarr store keeps counts in large compressed blocks. The Rust writer
targets about 4 million values per block, covering many groups. Each request
unpacks the whole block, keeps the one profile it asked for, and throws the
rest away. Asking for 300 groups x 20 length bins unpacks the same few blocks
6000 times.

Benchmark:

- Store: 300 groups x 20 length bins x 2001 positions, random `float32`,
  blocks of 99 x 20 x 2001 (the writer's 4M-value target), Zstd compression.
- Measured, 10 groups, first read in a fresh process: 3.5 s with the current
  code, 0.02 s with one request for all of them.
- *Inference*, all 300 groups: about 100 s with the current code, scaled from
  the 10-group run.

For a user, this is the "plot every TSS group" notebook cell taking minutes.

Fix:

Ask Zarr once for all selected groups and length bins together, so each block
is unpacked once. In zarr-python this is
`counts.oindex[group_indices, length_bin_indices, :]`. A probe confirmed it
returns unsorted selections such as groups `[5, 2]` and length bins `[3, 0]`
in the requested order. The public API does not change.

Tests: scalar, reversed, noncontiguous, and empty selections, keeping all
three dimensions and the requested order. The existing fixtures with
distinguishable values per axis are suitable. No timing thresholds in unit
tests.

### P-07 - Medium - Sparse data frames scan the whole file on every call

Evidence:

`py-cfdnalab/src/cfdnalab/ends.py:762-781` loops in Python over every stored
entry to find the requested rows and motifs.
`py-cfdnalab/src/cfdnalab/ref_kmers.py:754-773` is a copy of it.

Measured: `data_frame(window_idxs=0)` took 0.56 s on a store with 5.1M stored
entries. The cost grows with the size of the file, not the selection.

Impact:

Looping over windows costs the number of windows times the whole file.
Replacing the Python loop with a vectorized mask lowers the cost per call but
still scans every entry on every call. The R loader already does that, at 0.02
s per call on 1M entries. *Inference*: a loop over 50k windows would still take
about 17 minutes there.

Fix:

- Replace the Python loop with array lookups. This is the simple part.
- For repeated per-window access, find each row's range of entries with
  `np.searchsorted` on the row coordinates. This needs entries validated as
  sorted by row (ref-kmers validate this, end motifs do not, see P-29).
  Alternatively, keep a CSR index, at the cost of memory.
- Preserve the selector order after narrowing the entries.

### P-08 - High - Reference correction is slow and memory-hungry on genome-wide windows

Evidence:

`py-cfdnalab/src/cfdnalab/reference_correction.py:284-396` always builds long
pandas tables and joins them on text keys. The array and matrix outputs go
through the same path (`reference_correction.py:399-550`).

Measured, `sparse_corrected_counts_matrix(two_sided_correction="joint")`, 256
motifs, full support:

| Windows | Time | Peak after loading | Peak with correction |
| --- | --- | --- | --- |
| 3k | 2.2 s | not measured | not measured |
| 20k | 17.8 s | 0.36 GB | 2.32 GB |

*Inference*: genome-wide 10 kb windows (about 290k rows, the window size the
README recommends) would take about 4 to 5 minutes and need far more memory
than a laptop has.

Fix:

Compute the correction with matrices, like the Rust loader, and build the data
frame output from the result. The rewrite must keep these rules:

- Count positive support on the complete reference row, before matching it to
  sample columns. Reference motifs absent from the sample still count toward
  support and side frequencies. Example: frequencies `AA=1/4`, `CC=1/4`,
  `GG=1/2` and a sample count `_AA=6` give `6 / ((1/4) * 3) = 8`. Restricting
  the reference to the sample's columns first would give 24.
- Keep sparse references sparse and process bounded batches of rows. A dense
  rows x all-k-mers matrix recreates the memory problem for large k.
- Keep keyed row matching, global broadcasting, blacklist filtering on sample
  rows only, side aggregation before division, and first-occurrence order of
  the side axis.
- Unsupported zero counts stay zero. Unsupported positive counts raise an error
  or become `NaN`. `drop` stays available only for data frames.
- Corrected frequencies still need every sample motif of the row before motif
  selection, overflow-safe normalization, and `NaN` for the whole row when a
  positive count is unknown (`reference_correction.py:1210-1253`).

The scope of the unsupported-count check differs between loaders. Python and R
load every motif of the selected rows and fail on any unsupported positive
count, including motifs outside the motif selection. The Rust loader corrects
only the selected motifs (`src/output_loaders/reference_correction.rs:292-296`,
`:1010-1013`). A faster count-only path has to choose one of these rules on
purpose. Formula parity with Rust does not cover the frequency and selection
rules above.

Acceptance tests: reference and sample axes that differ, empty rows, and
recorded benchmark construction.

### P-09 - Low - Smaller costs

- Dense `sparse_counts_matrix(window_idxs=...)` reads the whole dense matrix
  before selecting rows (`py-cfdnalab/src/cfdnalab/ends.py:343-345`).
- Loading sparse ref-kmers checks entry order in a Python loop
  (`py-cfdnalab/src/cfdnalab/ref_kmers.py:1547-1554`). 0.9 s per 5M entries.
- The matrix path rebuilds a COO matrix from the stored arrays and converts it
  to CSR on every call (`ends.py:339-380`, same in `ref_kmers.py`). Whether an
  index or cache is worth it can be decided separately from P-07. The data
  frame helpers work on the COO arrays directly.

## API Consistency

### P-10 - Medium - Index column names are inconsistent

Data frames use `motif_index`, `length_bin`, and `position`. Selectors use
`motif_idx()`, `motif_idxs`, and `length_bin_idxs`. The `group_idx` and
`window_idx` columns do use `_idx`. A user writing `meta["motif_idx"]` gets a
`KeyError`.

The R package already uses `motif_idx`, `length_bin_idx`, and `position_idx`
columns, so the two packages also disagree with each other.

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
broad classes (for example short and long fragments) used to pick a profile.
There, midpoint and width add little.

Why it may still be worth adding:

A midpoint run with many narrow length bins gives a fragment length x position
heatmap (a V-plot). That needs a y-coordinate per length bin, which is
`length_midpoint_bp`.

Decision needed: keep the difference, add the two columns to midpoints, or
mention the difference in the midpoints docstring.

### P-12 - Medium - Corrected arrays do not say which rows they hold

Evidence:

- The windowed and grouped `corrected_counts_array` and
  `sparse_corrected_counts_matrix` accept `max_blacklisted_fraction`
  (`py-cfdnalab/src/cfdnalab/ends.py:1569`, `1676`, and the grouped versions).
- Rows above the cutoff are dropped
  (`py-cfdnalab/src/cfdnalab/reference_correction.py:441-447`, `500-506`).
- The docstrings do not mention the option (P-16).
- Column identity is available through `corrected_motifs_metadata()`. Row
  identity is not. `window_metadata()` and `group_metadata()` take no
  arguments.

Impact:

The filter is useful and should stay. After filtering, a user cannot easily
tell which windows the rows are. Plotting corrected values along the genome
means re-implementing the filter (`blacklisted_fraction <= cutoff`, applied
after the window selection, in selector order).

Fix:

- Document the option and the dropped rows on the corrected methods.
- Let `window_metadata()` and `group_metadata()` take the same row selector and
  `max_blacklisted_fraction`, and return rows in the same order. The returned
  indices can then be passed to the raw array methods, so adding the cutoff to
  every raw getter is optional.
- Row metadata for matrices must include rows that have no stored entries. It
  must not be derived from the distinct rows of a sparse corrected data frame.

### P-13 - Open - Raw end-motif data frames have no frequency column

A reason for leaving out the column was mentioned during review discussion but
is not recorded in the specs, the CHANGELOG, or the git history of `ends.py`.
Relevant facts:

- Ref-kmer frames have `frequency`. Corrected frames have
  `corrected_frequency`.
- `corrected_frequency` is normalized over the full stored motif axis before
  the motif selector is applied
  (`py-cfdnalab/src/cfdnalab/reference_correction.py:375-382`), so selecting
  motifs does not change it. The same rule could serve raw counts.
- For sparse output, the full stored axis is every observed motif, so the row
  total is still the total over all counted ends.
- *Inference* on possible reasons: for motifs-file runs, a frequency is
  relative to the targeted motifs, not to all fragment ends. For motif-group
  axes it is relative to the groups. Either could be misread as a share of all
  ends.

### P-14 - Low - Unknown labels raise different error types

An unknown motif raises `KeyError`, but an unknown side label in side mode
raises `ValueError` (`py-cfdnalab/src/cfdnalab/reference_correction.py:774`).

### P-15 - Low - `row_mode()` values follow CLI flag names

`EndMotifCounts.row_mode()` and `RefKmerFrequencies.row_mode()` return
`"size"`, `"bed"`, or `"grouped_bed"`. These mirror `--by-size`, `--by-bed`,
and `--by-grouped-bed` (`src/commands/cli_common.rs:408-437`,
`src/commands/ends/zarr.rs:287-292`). The classes are named `Windowed...` and
`Grouped...`. Lengths has no `row_mode()`, and its TSV does not record size
versus BED windows.

This is probably fine, since `"size"` versus `"bed"` carries information that
`"windowed"` would lose. At most, the `row_mode()` docstring could say which
flag each value comes from and which class it gives.

## Documentation

### P-16 - Medium - Corrected-count methods and all ref-kmer methods lack Parameters and Returns

This is not docstring inheritance:

- Each of the six `corrected_counts_array` and `sparse_corrected_counts_matrix`
  methods defines its own docstring
  (`py-cfdnalab/src/cfdnalab/ends.py:1126`, `1233`, `1561`, `1669`, `2034`,
  `2143`). None has a Parameters or Returns section. `allow_densify`,
  `max_blacklisted_fraction`, `use_global_bias`, `unsupported_motifs`, and
  `two_sided_correction` are undocumented there.
- All 32 public ref-kmer methods define their own docstrings and none has
  Parameters or Returns, including those that take arguments (`data_frame`,
  `motif_idx`, `group_idx`, and the array and matrix methods). `read_ref_kmers`
  (`py-cfdnalab/src/cfdnalab/ref_kmers.py:1287`) has a 3-line docstring.
  `read_end_motifs`, `read_lengths`, and `read_midpoints` have both sections.
- The website generator parses each class body statically with
  `ast.get_docstring` (`website/scripts/generate_loader_docs.py:789-829`), so
  nothing is inherited there either. The website pages have the same gaps.

Suggestion: `validate_python_docstrings`
(`website/scripts/generate_loader_docs.py:1301-1327`) only checks that a
docstring exists. It could also require a Parameters section when the method
takes arguments.

### P-17 - Low - The reference-correction explanation is copied into 9 docstrings

The about 60-line "Reference correction" and "Two-sided correction modes"
sections appear in 9 docstrings in `ends.py`. The copies are identical today.

Ways to keep them from drifting:

- Keep the copies and add a test that compares those sections across the 9
  docstrings. No change for readers. Cheapest.
- Write the full explanation once (for example in the website ends guide) and
  replace each copy with a short summary and a pointer. Shorter docstrings, but
  `help(...)` no longer shows the full explanation.
- Build the docstrings at import time from one shared string. `help(...)` would
  work, but the website generator reads source statically and would show the
  unexpanded text. Not recommended without generator changes.

Recommended start: the test.

### P-18 - Notes - README and guides

- "was used" should be "were used". "why we recommend" should be "which is why
  we recommend".
- The eligible-intervals example (`py-cfdnalab/README.md:78-88`) loops group by
  group with `iterrows`. One call with `group_idxs=` does the same.
- `website/docs/guides/fragment_lengths_guide.md:119-145` calls
  `group_metadata()` and `window_metadata()`, and selects by `groups=` and
  `window_idxs=`, on the same object. A file supports only one of these. Label
  them as alternatives per file mode.
- End-to-end recipes would help: a TSS profile normalized to its flanks, a
  short/long ratio per window, a samples x motifs frequency matrix.

### P-19 - Low - Length counts are called "raw"

`counts_array` says "raw length counts". Stored counts already include the
run's window assignment and optional weighting:

- GC correction and scaling factors weight counts when enabled.
- Count-overlap window assignment gives fractional counts even without weights
  (`src/commands/lengths/lengths.rs:1247-1250`,
  `src/commands/ends/ends.rs:1385-1388`).

`length_settings.json` records `gc_correction_used` and `scaling_factors_used`
(`src/commands/lengths/writer.rs:75-76`). `end_settings.json` records neither
(`src/commands/ends/write.rs:53-87`), so a settings accessor would not tell
end-motif users whether counts were weighted.

Suggestion: describe them as stored counts with the run's window assignment and
weighting already applied.

## Tests

- P-20: The midpoint fixture uses positions the writer never produces (P-04).
- P-21: No test covers numeric-looking labels: all-numeric group names,
  zero-padded chromosomes, or `1.10` next to `1.1` (P-01).
- P-22: Corrected values for motif-group outputs are not tested anywhere: not
  in the Python unit tests, the downstream tests, or the Rust loader tests.
  Only the rejection of two-sided modes is tested. A useful case: group labels
  containing underscores (so labels must match verbatim), sample and reference
  columns in different orders, and reference support beyond the selected
  groups. With counts 6 and 4 and group frequencies 1/4 and 3/4, the corrected
  counts are 12 and 8/3 and the corrected frequencies are 9/11 and 2/11. The
  docs should say that this normalizes over supported groups, whose membership
  must match between runs.
- P-23: Data frames already have tests for blacklist filtering and `NaN`
  propagation. The gap is at the fixed-shape array and matrix boundary:
  - a reordered row selection followed by blacklist filtering,
  - a retained all-zero sparse row,
  - explicit `NaN` storage in the sparse matrix under `keep_na`.

  Check values and axis identity together.
- P-24: Unit fixtures mostly use end-motif schema version 1 without
  `motif_axis_kind`. The writer produces version 2. The downstream tests cover
  version 2.
- `tests/test_package.py` only checks that public symbols exist.
- The midpoint group name `"gamma-unicode-aa"` suggests a Unicode test, but the
  name is ASCII.

## Design Questions

### P-25 - Low - Empty sparse end-motif stores cannot be opened

`py-cfdnalab/src/cfdnalab/ends.py:2308-2327` rejects a sparse store with no
stored counts, so even its window metadata is unreadable. A dense store with
zero counts opens fine.

This is a deliberate choice: CHANGELOG 0.3.0 records the explicit no-count
error, and tests assert it. The R loader has to reject these stores, because
the CRAN `zarr` reader builds the whole hierarchy when opening and cannot handle
zero-length arrays (`r-cfdnalab/R/end_motifs.R:226-235`).

Open question: whether Python should open such stores for metadata only. A
Python-only change has to define selection and correction on an empty axis,
where side widths cannot be inferred from labels.

### P-26 - Low - Blacklist cutoffs can keep different windows in lengths and ends

Lengths TSVs round `blacklisted_fraction` to 3 decimals
(`src/commands/lengths/writer.rs:24`). Zarr outputs keep full precision. A
window with `0.1004` blacklisted is `0.100` in lengths and passes a `0.1`
cutoff there, but fails it in ends. This assumes the same windows and
blacklist. Grouped fractions are only comparable when the eligible interval
sets also match.

Fix in the writer (keep more precision) and document the difference. A loader
tolerance would change the requested cutoff and cannot recover the discarded
digits.

### P-27 - Medium - Complement-collapsed output gets a wrong reference correction

Evidence:

- `--collapse-complement` is accepted in ordinary builds. `hide` only removes
  it from the help text (`src/commands/ends/config.rs:225-242`).
- The settings entry is only written when the `ends_experimental` feature is
  enabled (`src/commands/ends/write.rs:126-137`). Neither the Zarr root nor,
  in ordinary builds, `end_settings.json` records the collapse.
- Collapsing merges each motif with its same-orientation complement, not its
  reverse complement (`src/commands/ends/output.rs:86-93`).

Impact:

With `_A` and `_T` merged into one label, the matching reference mass is
`p(A) + p(T)`, with support counted on the collapsed axis. The current
correction uses `p(A)` and uncollapsed support. The result is wrong, and no
error is raised. Exposure is low (hidden flag, documented as possibly not
biologically meaningful), but the guard is cheap, so it belongs before
validation.

Fix:

Record the collapse in the Zarr root. Make Python, R, and Rust correction
reject collapsed stores until they implement the same transformation on the
reference side.

### P-28 - Low - Some errors are hidden or misreported

- `_has_array` swallows every exception, so a damaged store is reported as a
  missing array. It exists in `ends.py:2446-2454`, `midpoints.py:580`, and
  `ref_kmers.py:1483`.
- The `zstd` command-line fallback in `_read_length_counts_table`
  (`py-cfdnalab/src/cfdnalab/lengths.py:784-822`) runs after any pandas error
  on a `.zst` file, including parse errors. It cannot help when `zstandard` is
  missing, because the module imports that required dependency and the header
  reader uses it first. It mostly replaces the original pandas error with a
  generic message. Remove it, or narrow it, and keep the parse or
  decompression error in the message.

### P-29 - Low - Validation gaps between the readers

These can be fixed with targeted checks. The shared selector and blacklist
helpers in `_helpers.py` already exist, so a new shared module is not needed.

- End-motif sparse arrays are not checked for sorted, unique coordinates.
  Ref-kmers check this. Duplicate coordinates show up as separate data frame
  rows but are summed when converted to a sparse matrix.
- End-motif window intervals and blacklist fractions are not validated
  (`ends.py:199-221`, `:2469-2525`). Ref-kmers validate both.
- Midpoints validate shapes and axes (`midpoints.py:137-152`) but not
  coordinate integer types or valid bin intervals. The negative-position
  fixture (P-20) loads because of this. Finite negative count values must stay
  valid, because smoothing produces them.
- The R end-motif sparse reader has the same missing order check (R-11).

Specify each malformed-store rejection as a test. Whether to extract shared
code can be decided later.
