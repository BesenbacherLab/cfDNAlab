# r-cfdnalab Loader Review

Date: 2026-11-05

Scope: `r-cfdnalab/R/*`, `r-cfdnalab/tests/testthat/*`, `README.md`,
`DESCRIPTION`, `NAMESPACE`, `DEVELOPMENT.md`, `NEWS.md`, the reference
correction tests in `downstream_tests/R/test_cfdnalab_r_package.R`, and the R
snippets in `website/docs/guides/`.

Method: all R sources were read, including the relevant internals of the CRAN
`zarr` reader. Claims about runtime behavior, performance, and locale handling
were checked with probes using R 4.4.1, CRAN `zarr` 0.4.0, data.table 1.18.4,
and `pkgload::load_all()`. Probe stores used the Rust writer's fill values and
chunk layout. Memory figures are process peak resident memory from
`/usr/bin/time -l`, reported next to the peak after loading only. The test
suite was not run. Static documentation checks (`tools::checkDocFiles`,
`tools::codoc`, `tools::undoc`) were run, and `roxygen2::roxygenise()` was run
in a scratch copy to compare the generated docs. Anything marked *inference*
was not measured.

Findings shared with the Python package are not repeated here. R-11 points to
the matching entries in `python_loader_review.md`.

## Checked And Not Findings

- R and Python give the same corrected values (to 10 decimals) for `joint`,
  `split`, `outside`, and `inside` on the same 3k-window store.
- Every hand-derived expectation in
  `tests/testthat/test-reference-correction-core.R` and in the two two-sided
  correction tests in `downstream_tests/R/test_cfdnalab_r_package.R:116-334` is
  correct.
- The one-based indexing rules in `DEVELOPMENT.md` are followed. Public data
  frames use one-based `*_idx` columns and never expose `*_idx0`.
- Static documentation checks report no undocumented arguments and no usage
  mismatches.
- Midpoint reads are fast enough. The `zarr` reader keeps decoded blocks, so
  each extra profile costs about 1 ms. `midpoint_data_frame()` for 150 groups x
  20 length bins (6M rows) took 4.5 s. Python P-06 does not apply.
- Sparse stored-row data frames use vectorized lookups
  (`R/end_motifs.R:1377-1390`, `R/ref_kmers.R:1376-1389`): 0.02 s for one
  window on a 1M-entry store. Each call still scans every stored entry, so the
  repeated-access part of Python P-07 applies.
- The sparse ref-kmer order check is vectorized (`R/ref_kmers.R:557-573`).
- Other `order()` calls sort integers. `sort()` on labels in
  `R/reference_correction.R:1237` and `1263` only builds error messages.
- Data frame index columns already use `*_idx` names consistently
  (`motif_idx`, `length_bin_idx`, `position_idx`).

## Backend Behavior That Constrains Fixes

Two properties of `zarr_array$read()` in CRAN `zarr` 0.4.0 shape several
findings below.

- Selection: `read(selection = list(...))` takes the minimum and maximum of each
  axis selection and reads that whole range in file order. Gaps and the
  requested order are ignored. Probe: selecting groups `c(3, 1)` returns groups
  1, 2, and 3. Selective reads must read a bounded range and then pick the
  requested indices in order, or make several bounded reads.
- Fill values: integer arrays mark exact matches to the fill value as `NA`.
  Floating arrays mark every value within a tolerance of the fill value as
  `NA`. The tolerance is `max(1.49e-08 * max(abs(slice), abs(fill)), 1e-12)`,
  where `max(abs(slice))` is the largest absolute value in the whole returned
  slice. A `NaN` fill marks only `NaN` values. There is no option to turn this
  off.

## Release Triage

Fix before validation:

- R-08: Valid floating values can become `NA`, and some stores fail to load.
- R-01: Canonical ref-kmer outputs fail to load under Danish collation.
- R-02: `end_motif_data_frame(densify = TRUE)` builds the dense matrix for the
  whole store.
- R-03: Numeric-looking chromosome and group labels are changed or rejected.
- R-04: Reading `.tsv.zst` needs a separately installed `zstd` program.
- R-05: README and website examples that fail or name the wrong flag.

Fix before people use it on genome-wide data:

- R-06: Reference correction is slow and memory-hungry.

Docs and API:

- R-07: Count matrices have no row or column names.
- R-09: Committed man pages are out of date.
- R-10: No way to get a selected 3D midpoint array.

## Wrong Results Or Broken Use

### R-08 - High - Valid floating values can become NA, and some stores fail to load

Evidence:

- The fill-value rule above applies to every floating array the R loader
  reads. The Rust writer uses `-1` as the fill value for floating arrays
  (`src/shared/zarr.rs:49-55`) and chooses fill values outside each array's
  data domain so that readers like this one do not mask real data
  (`.AI/docs/specs/zarr_schema.md:21-29`).
- Midpoint counts are smoothed by default and can be negative (P-04), so `-1`
  is inside their data domain. `src/shared/zarr.rs:49-51` describes them as
  non-negative.
- Because the tolerance grows with the largest value in the slice, non-negative
  arrays are also affected once that value is large enough. At about `6.7e7`
  the tolerance reaches 1, and zero then falls inside it.

Probes:

- `[2^27, 0, 1, 0.5, 3.5]` with fill `-1` reads back as
  `[2^27, NA, NA, NA, 3.5]`.
- A midpoint-like array where one group has a peak of `1e7` and another group
  has a smoothed value of `-1.1`: reading the second group alone returns
  `-1.1`, and reading the whole array returns `NA` for that value. The same
  value depends on which slice is read.
- A sparse end-motif store with one stored count of `2^27` and other counts of
  1 to 49 fails `read_end_motifs()` with
  `sparse/count must contain finite non-negative values`, because the small
  counts were turned into `NA`.
- Exact `-1` values are reachable from integer counts. With
  `--smoothing savgol=5`, the neighborhood `[23, 0, 2, 0, 0]` gives
  `(-3*23 + 17*2) / 35 = -1` (`src/commands/midpoints/smoothing.rs:107-116`).

Impact:

- Midpoint values near `-1` become `NA` without a warning. How near depends on
  the largest value in the slice that is read.
- Count, frequency, and scaling-factor arrays lose small values when the same
  array also holds a value around `6.7e7` or more. For sparse arrays read at
  load time, this makes the load fail. For dense arrays read later, it gives
  `NA`. *Inference*: mostly global outputs from deep or pooled data with
  weighted fractional counts. Rare otherwise.
- Python is not affected.

Fix (Rust side):

Use `NaN` as the fill value for floating arrays and update the schema spec and
the writer comments. Do not use `-Inf`. With an infinite fill, the tolerance
becomes infinite and every finite value is masked. Probe: `[0, 1, -0.5, 2]`
with fill `-Inf` reads back as all `NA`. A `NaN` fill masks only `NaN`, and
`[0, -1, -0.5, 2]` with fill `NaN` reads back unchanged.

Outputs already written with a finite fill cannot be repaired by the R loader,
because the values are lost inside `read()`.

Tests: writer-backed downstream fixtures whose metadata has a `NaN` fill,
asserting that finite negative midpoint values (including `-1`) survive full
and sliced reads, including a shortened last chunk.

### R-01 - High - Canonical ref-kmer outputs fail to load under Danish collation

Evidence:

- `R/ref_kmers.R:515-526` picks the canonical even-k label with
  `min(motif, reverse_complement)`. R compares strings with the current
  collation locale, not by bytes.
- Danish collation treats the pair `"AA"` as a single letter (A-ring), which
  sorts after `"Z"` and therefore after `C`, `G`, and `T`. Labels never contain
  `Z`, but that does not matter. Any comparison whose first difference lands on
  an `"AA"` goes the wrong way. Probe: `min("AACG", "CGTT")` is `"AACG"` under
  `C` and `"CGTT"` under `da_DK.UTF-8`, because `"AACG"` compares as A-ring
  followed by `"CG"`.
- End to end: a canonical store with `AACG`, `ACGT`, and `CCGG` loads under
  `C` and fails under `da_DK.UTF-8` with
  `canonical reference k-mer motif label 'AACG' at motif_idx 1 should be 'CGTT'`.
- testthat edition 3 forces `C` collation during tests
  (`testthat::local_reproducible_output()` calls `withr::local_collate("C")`),
  so the test suite cannot catch this.

Impact:

Canonical labels that R computes differently from Rust under `da_DK.UTF-8`,
counted over full k-mer sets:

| k | Canonical labels | Wrong under Danish collation | Example |
| --- | --- | --- | --- |
| 2 | 10 | 1 | `AA` (versus `TT`) |
| 4 | 136 | 18 | `AACG` (versus `CGTT`) |
| 6 | 2,080 | 298 | |
| 8 | 32,896 | 4,806 | |

About 15% of labels are affected. A load fails when the stored axis contains at
least one affected label. A complete even-k axis always does. A sparse axis
with only unaffected labels loads. Any R session with Danish collation is
affected. The macOS system locale on the development machine is `da_DK`, so
sessions that take their collation from the system locale are affected
(*inference*: R GUIs and RStudio on macOS do). Norwegian collation has the same
`"AA"` rule (*inference*, not probed). Odd k uses the middle-base rule and is
not affected.

This is the only place in the R package where locale ordering changes a
result. Python compares strings by code point and is not affected.

Fix:

Compare in byte order with radix ordering, which sorts strings in byte order
whatever the locale is. This is correct for the validated `ACGT` alphabet. In
`cf_canonical_ref_kmer()` (`R/ref_kmers.R:517-519`):

```r
if (nchar(motif, type = "bytes") %% 2L == 0L) {
  # Compare in byte order like the Rust writer. Locale collation can reorder
  # labels, e.g. Danish collation treats "AA" as one letter after "Z"
  candidates <- c(motif, reverse_complement)
  return(candidates[order(candidates, method = "radix")[1L]])
}
```

Verified under `da_DK.UTF-8`: for all 4,096 6-mers this picks exactly the
labels Rust picks. `sort(c("CGTT", "AACG", "TAAA"), method = "radix")` gives
`AACG CGTT TAAA`, while the default `sort()` gives `CGTT TAAA AACG`.

Setting the collation instead is also correct. `withr::local_collate("C")` at
the top of a function sets `C` collation and restores the session's setting
when the function exits, including on errors (checked with withr 3.0.2).
`withr` is already installed with testthat and would become an `Imports`
dependency. Radix ordering is the smaller change for this single comparison.
`withr::local_collate("C")` is the better convention if more functions start
ordering or comparing labels.

Tests:

Add `withr` to `Suggests` and set Danish collation inside the test with
`withr::local_collate("da_DK.UTF-8")`. It is restored automatically when the
test ends, also when an assertion fails. Then assert fixed expected labels
(`AA` versus `TT`, `AACG` versus `CGTT`) and load a canonical fixture containing
`AACG`. Both fail with the current code.

`withr::local_collate()` only warns when a locale is not installed and keeps
the previous collation. Probe: requesting `xx_XX.UTF-8` left `en_US.UTF-8`
active and `min("AACG", "CGTT")` returned `"AACG"`. The test must therefore
check that Danish collation is really active, for example by skipping unless
`min("AACG", "CGTT") == "CGTT"` after setting it. Otherwise it passes without
testing anything on machines without the locale. Locale names differ by
platform, so try `da_DK.UTF-8` and the Windows name `Danish_Denmark.1252`
before skipping.

### R-02 - High - `end_motif_data_frame(densify = TRUE)` densifies the whole store

Evidence:

- `R/end_motifs.R:1341-1345` calls `dense_counts_matrix(x, allow_densify =
  densify)` without the selected rows and motifs, then subsets the result.
- Reference correction goes through the same function. That covers corrected
  data frames with `densify = TRUE` and every `dense_corrected_counts_matrix()`
  call (`R/reference_correction.R:277-283`, `436`).
- The ref-kmer version passes the selection (`R/ref_kmers.R:1333-1338`), so
  this is a bug and not a design choice.

Measured on a sparse store with 50k windows x 1024 motifs (1M stored entries):

- `end_motif_data_frame(window_idxs = 1, motifs = <2 motifs>, densify = TRUE)`:
  about 405 MB extra peak memory. The full dense matrix is 391 MB.
- The same selection through `dense_counts_matrix(..., window_idxs =,
  motifs =)`: 45 MB.

Impact:

The README example at `README.md:112-118` does exactly this. With genome-wide
windows and larger motif sets, a request for a handful of rows runs out of
memory. *Inference*: 290k windows x 1024 motifs is 2.4 GB. A 4+4 motif axis is
far beyond any laptop.

Fix:

Two separate changes:

1. Sparse stores: call `cf_dense_end_motif_matrix_for_indices(x, row_indices,
   motif_indices, densify)` so only the selection is densified.
2. Dense stores: the dense branch of that helper still reads the whole counts
   array before subsetting (`R/end_motifs.R:632-638`, and `533` for
   `sparse_counts_matrix()` on dense stores). Ref-kmers do the same in
   `cf_ref_kmer_dense_frequency_matrix_for_indices()`
   (`R/ref_kmers.R:986-994`). Making these reads selective must follow the
   backend selection rule above: read a bounded range, then pick the requested
   indices in order, keeping length-one dimensions.

Reference correction still loads every motif of the selected rows on purpose,
because corrected frequencies need the whole row (P-08). Fix 1 bounds the rows,
not the motif axis.

### R-03 - High - Numeric-looking chromosome and group labels are changed or rejected

Same bug as Python P-01.

Evidence:

- `R/lengths.R:454` and `457` call `data.table::fread()` without column
  types, so it guesses.
- `R/lengths.R:459-461` converts `chrom` back to text after parsing. By then
  the original text is gone.
- `R/lengths.R:595` rejects a `group_name` column that is not character.

Probes:

- Windowed chromosomes `001` and `1` both load as `"1"`. No error.
- Group names `1`, `2` and `001`, `002` and `1.10`, `1.1` all fail with
  `group_name must contain non-empty character strings`.

Fix:

Read `chrom` and `group_name` as character at parse time with `colClasses`.
Pass only the columns present in the file, because `fread()` warns about absent
names. Remove the cast after parsing. Use the same test cases as P-01, for
plain and `.tsv.zst` input.

### R-04 - Medium - Reading `.tsv.zst` needs a separately installed `zstd` program

Evidence:

- `R/lengths.R:446-454` calls the `zstd` command-line program and stops if it
  is missing. `.tsv.zst` is the default `cfdna lengths` output.
- The requirement is not mentioned in `README.md`, `DESCRIPTION` (no
  `SystemRequirements`), or `website/docs/guides/fragment_lengths_guide.md`.
- `qs2` is already a required dependency and exports `zstd_decompress_file()`.
  Probe: it decodes a zstd file written in streaming mode (no decompressed size
  in the header, like the Rust writer's `zstd::Encoder` at
  `src/shared/writers.rs:21`). `qs2::zstd_decompress_raw()` fails on the same
  file with `negative length vectors are not allowed`, so use the file version.

Impact:

*Inference*: macOS and Windows do not ship `zstd`, so many R users cannot read
lengths output without installing an extra program first.

Fix:

Decompress with `qs2::zstd_decompress_file()` to a temporary file, then
`fread()` it. Drop the command-line path. Requirements:

- Works on streaming frames without a declared size, tested against real Rust
  writer output.
- Reports corrupted and truncated input clearly.
- Removes the temporary file on success and on error.
- Declares a minimum `qs2` version. The installed 0.1.7 has
  `zstd_decompress_file()`. When it was added is not verified (the installed
  package has no NEWS file).

The compressed-read part of the current test only runs when `zstd` is installed
(`tests/testthat/test-lengths.R:511-527`). The fixture construction must stop
needing the program, so compressed loading is tested everywhere.

### R-05 - Medium - README and website examples that fail or name the wrong flag

- `README.md:235-239` and `website/docs/guides/fragment_ends_guide.md:301-305`
  call `dense_corrected_counts_matrix(...)` without `allow_densify = TRUE`.
  This errors on default sparse `cfdna ends` output
  (`R/reference_correction.R:420-426`). The guide's `--all-motifs` is passed to
  `cfdna ref-kmers` and does not make the sample output dense. Same as Python
  P-02. Prefer `sparse_corrected_counts_matrix()` in general examples.
- `README.md:80` says `"size"` rows come from `--window-size`. The flag is
  `--by-size` (`src/commands/cli_common.rs:408`).
- `website/docs/guides/fragment_lengths_guide.md:158-163` calls
  `window_metadata(lengths)` and `group_metadata(lengths)` one after the
  other. Only one exists for a given file, so the block stops with an error.
  The Python block at `:119-145` has the same problem (P-18). Label these as
  alternatives per file mode.
- Notes only: `README.md:207` "why we recommend" should be "which is why we
  recommend". `README.md:209` "was used" should be "were used".

## Performance At Realistic Sizes

### R-06 - High - Reference correction is slow and memory-hungry

Evidence:

`R/reference_correction.R` builds long data frames and uses:

- base `merge()` (`945`, `1069`, `1197`),
- `stats::aggregate()` (`1064`, `1161`),
- `stats::ave()` with an R function per output row (`1318`),
- an R function call per row to split labels (`cf_add_end_motif_sides()` at
  `1111-1114` calls `cf_split_end_motif_label()`, which runs `gregexpr()`, on
  every row).

Measured, `sparse_corrected_counts_matrix()`, 256 motifs, full support:

| Windows | Mode | R time | Python time | R peak (load only / with correction) | Python peak (load only / with correction) |
| --- | --- | --- | --- | --- | --- |
| 3k | joint | 6.7 s | 2.2 s | not measured | not measured |
| 3k | split | 12.4 s | not measured | not measured | not measured |
| 3k | outside | 10.4 s | not measured | not measured | not measured |
| 20k | joint | 70.9 s | 17.8 s | 1.12 GB / 3.39 GB | 0.36 GB / 2.32 GB |

R takes about 4 times as long. The extra memory for the correction itself is
similar (about 2.3 GB in R and 2.0 GB in Python). R's higher total comes from
its larger baseline after loading.

*Inference*: genome-wide 10 kb windows (about 290k rows) would take about
15 to 20 minutes for `joint` and need far more memory than a laptop has.

Fix:

Same approach and the same rules as Python P-08, using `Matrix` and R's
one-based public indices:

- Count positive support on the complete reference row before matching it to
  sample columns.
- Keep sparse references sparse and process bounded batches of rows. Avoid a
  dense rows x all-k-mers intermediate.
- Split motif labels into sides once per unique label, not once per row.
- Keep the unsupported-count, row matching, side-axis order, and frequency
  rules listed in P-08.

Do R-02 separately and first, so the memory effect of each change can be
measured on its own.

## API And Documentation

### R-07 - Medium - Count matrices have no row or column names

Evidence:

- `sparse_counts_matrix()` and `dense_counts_matrix()` for end motifs and
  ref-kmers return matrices without dimnames (`R/end_motifs.R:530-549`,
  `626-646`, `R/ref_kmers.R:780-797`, `858-881`). Probe on the downstream
  sparse windowed fixture confirms `NULL` dimnames.
- Corrected matrices get column names only in `outside` and `inside` modes
  (`R/reference_correction.R:446-448`, `547-551`). Unlike Python, R has no
  `corrected_motifs_metadata()`, so in `joint` and `split` mode the columns
  are only identified by matching the selector order.
- Length count matrices have `count_*` column names (`R/lengths.R:639`).
- `dense_counts_vector()` and `dense_frequencies_vector()` are named by motif
  (`R/end_motifs.R:975`, `R/ref_kmers.R:1094`, `1106`).

Impact:

After `motifs = c("AC_GT", ...)` the user gets a matrix that does not say which
column is which motif. Names would not align samples automatically: `rbind()`
and `cbind()` on matrices bind by position. They do make the columns
inspectable and allow explicit reindexing (`m[, shared_motifs]`) before
binding, which is how misalignment from P-03 is avoided.

Fix:

- Name the columns of raw and corrected matrices in every mode, including empty
  results.
- Keep the public `window_idx` or `group_idx` in row identity, for example as
  row names or in returned metadata. Coordinates alone may not identify a BED
  row uniquely. No duplicate-interval check was found in
  `src/shared/bed.rs`, and the R correction code already includes
  `window_idx` in its row key. This matters most after blacklist filtering
  (P-12).
- Optionally add dimnames to `midpoint_array()` and `profile_array()` (group
  names, length bins, position starts).
- Document reindexing to a shared motif order before binding samples.

Tests: the existing `unname()` value assertions
(`downstream_tests/R/test_cfdnalab_r_package.R:198`, `217`, `298`, `327`) can
stay. Add separate assertions for the exact dimnames and their order.

### R-09 - Low - Committed man pages are out of date

Evidence:

Running `roxygen2::roxygenise()` (8.0.0, the version in `DESCRIPTION`) in a
scratch copy changes three files:

- `man/dense_corrected_counts_matrix.Rd`
- `man/end_motif_data_frame.Rd`
- `man/sparse_corrected_counts_matrix.Rd`

Each gains the paragraph saying that correction requires `--orientation both`,
is approximate, and needs windows of at least a few kilobases.
`website/docs/generated/loaders/r-api.mdx` is generated from the Rd files and
does not mention `--orientation both`. `NAMESPACE` is in sync.

Fix:

Regenerate the man pages and the website pages. Consider a check that fails
when `roxygenise()` changes tracked files.

### R-10 - Medium - No way to get a selected 3D midpoint array

Evidence:

- `midpoint_array()` reads the whole store (`R/midpoints.R:395-398`).
- `profile_array()` returns one profile (`R/midpoints.R:190-210`). Its
  singular selector names fit a single-profile function.
- Python's `counts_array()` takes the same group and length-bin selectors as
  its data frame method.

Impact:

A heatmap of 50 selected groups means looping over `profile_array()` or
reading everything. For a large grouped store, reading everything is not an
option (*inference* on sizes).

Fix:

Give `midpoint_array()` the selectors from `midpoint_data_frame()`. Because of
the backend selection rule, a single read for groups `c(5, 2)` returns groups 2
to 5 in file order. Read the bounded range, then pick the requested groups and
length bins in order. A contiguous selection of 10 groups x 20 length bins x
2001 positions read in one call took 0.01 s.

Tests: noncontiguous and reversed group and length-bin selections, scalar
axes, and empty selections that make no backend read. Assert values with the
existing fixtures that have distinguishable values per axis, not just shapes.

### R-11 - Shared with the Python review

These apply to R as well. Details are in `python_loader_review.md`.

- P-03: Motif axes differ between samples. `motifs()` and `motif_idx` are
  local to each file.
- P-04: `read_midpoints()` (`R/midpoints.R:19-35`) does not explain values,
  smoothing, normalization, or positions either.
- P-05: R docs explain motif direction only in the correction context
  (`R/end_motifs.R:985-986`).
- P-11: The `length_bins()` column difference exists in R too
  (`R/lengths.R:67-77` versus `R/midpoints.R:131-139`).
- P-12: Corrected matrices filtered by `max_blacklisted_fraction` drop rows
  (`R/reference_correction.R:427-428`, `495-496`). With no row names (R-07),
  nothing says which rows remain.
- P-13: R end-motif data frames have no frequency column either.
- P-19: `R/lengths.R:252` also calls the counts "raw".
- P-22: Motif-group correction values are not tested in R either.
  `tests/testthat/test-reference-correction.R:124-141` only checks that
  two-sided modes are rejected.
- P-25: R has to reject empty sparse stores, because the `zarr` reader cannot
  open zero-length arrays (`R/end_motifs.R:226-235`). Supporting them would
  mean bypassing the reader's eager hierarchy construction.
- P-26: Blacklist rounding in lengths TSVs.
- P-27: Complement collapse is not recorded, and R correction would use the
  wrong reference mass as well.
- P-29: End-motif sparse coordinates are not checked for sorted, unique order
  (`R/end_motifs.R:136-147`), while ref-kmers check it
  (`cf_validate_ref_kmer_sparse_coordinates()`). This needs no refactor.

### R-12 - Low - Loading large k-mer axes is slow

Ref-kmer label validation runs an R loop per label (`R/ref_kmers.R:447-492`).
Motif ASCII decoding calls `apply()` per row (`R/helpers.R:1031-1033`).

Measured `read_ref_kmers()` on global stores with every k-mer:

| k | k-mers | R | Python |
| --- | --- | --- | --- |
| 8 | 65,536 | 0.8 s | 0.1 s |
| 10 | 1,048,576 | 9.9 s | 1.2 s |

Fine for k up to 8.

### R-13 - Low - Package metadata

- `DESCRIPTION` sets `LazyData: true` without a `data/` directory. *Inference*:
  `R CMD check` reports this as a NOTE (not run).
- No `SystemRequirements: zstd`. This goes away if R-04 is fixed.
- No minimum `qs2` version (R-04).

## Tests

- No test runs under a non-`C` collation (R-01). testthat forces `C`
  collation, so the test must set Danish collation itself with
  `withr::local_collate()` and skip when it does not take effect.
- No test covers numeric-looking labels (R-03).
- No test would catch R-02, because the result is correct. Do not rely on a
  large declared shape and allocation failure, which may succeed or may kill
  the R process. Keep fixtures tiny, and replace the whole-store step with a
  version that fails when called, for example with
  `testthat::local_mocked_bindings()` on `dense_counts_matrix`. Then a bounded
  selection with `densify = TRUE` must still return the expected values, order,
  and retained empty rows. Test the dense selective-read fix separately.
- The compressed-read part of the lengths test only runs when `zstd` is
  installed (R-04).
- Corrected matrix comparisons strip names with `unname()`. Exact dimnames are
  not asserted anywhere (R-07).
- No test covers R-08. It needs writer-backed fixtures.
- Motif-group correction values are untested (R-11).
