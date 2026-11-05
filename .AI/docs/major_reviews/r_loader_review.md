# r-cfdnalab Loader Review

Date: 2026-11-05

Scope: `r-cfdnalab/R/*`, `r-cfdnalab/tests/testthat/*`, `README.md`,
`DESCRIPTION`, `NAMESPACE`, `DEVELOPMENT.md`, `NEWS.md`, the reference
correction tests in `downstream_tests/R/test_cfdnalab_r_package.R`, and the R
snippets in `website/docs/guides/`.

Method: all R sources were read. Claims about runtime behavior, performance,
and locale handling were checked with probes using R 4.4.1, CRAN `zarr` 0.4.0,
data.table 1.18.4, and `pkgload::load_all()`. Probe stores were written with
the same fill values and chunk layout as the Rust writer. The test suite was
not run. Static documentation checks (`tools::checkDocFiles`, `tools::codoc`,
`tools::undoc`) were run, and `roxygen2::roxygenise()` was run in a scratch
copy to compare the generated docs. Anything marked *inference* was not
measured.

Findings shared with the Python package are not repeated here. R-11 points to
the matching entries in `python_loader_review.md`.

## Checked And Not Findings

- R and Python give the same corrected values (to 10 decimals) for `joint`,
  `split`, `outside`, and `inside` on the same 3k-window store.
- Every hand-derived expectation in `tests/testthat/test-reference-correction-core.R`
  and in the two two-sided correction tests in
  `downstream_tests/R/test_cfdnalab_r_package.R:116-334` is correct.
- The one-based indexing rules in `DEVELOPMENT.md` are followed. Public data
  frames use one-based `*_idx` columns and never expose `*_idx0`.
- Static documentation checks report no undocumented arguments and no
  usage mismatches.
- Sparse stored-row data frames are built with vectorized lookups
  (`R/end_motifs.R:1377-1390`, `R/ref_kmers.R:1376-1389`). 0.02 s for one
  window on a 1M-entry store. Python P-07 does not apply.
- Midpoint reads are fast enough. The `zarr` reader keeps decoded blocks, so
  each extra profile costs about 1 ms. `midpoint_data_frame()` for 150 groups x
  20 length bins (6M rows) took 4.5 s. Python P-06 does not apply.
- The sparse ref-kmer order check is vectorized (`R/ref_kmers.R:557-573`).
- Other `order()` calls sort integers. `sort()` on labels in
  `R/reference_correction.R:1237` and `1263` only builds error messages.
- Data frame index columns already use `*_idx` names consistently
  (`motif_idx`, `length_bin_idx`, `position_idx`).

## Release Triage

Fix before validation:

- R-08: Valid midpoint values can silently become `NA`. The annotation below
  also identifies an unsafe proposed replacement fill value.
- R-01: Canonical ref-kmer outputs fail to load under Danish collation.
- R-02: `end_motif_data_frame(densify = TRUE)` builds the dense matrix for the
  whole store.
- R-03: `read_lengths()` fails when every group name is numeric.
- R-04: Reading `.tsv.zst` needs a separately installed `zstd` program.
- R-05: README and website examples that fail or name the wrong flag.

Fix before people use it on genome-wide data:

- R-06: Reference correction is slow and memory-hungry.

Docs and API:

- R-07: Count matrices have no row or column names.
- R-09: Committed man pages are out of date.
- R-10: No way to get a selected 3D midpoint array.

## Wrong Results Or Broken Use

### R-01 - High - Canonical ref-kmer outputs fail to load under Danish collation

> **Review note - scope and portable coverage:** Failure requires an affected
> stored label. A complete canonical even-k axis contains one, but a sparse
> axis containing only unaffected labels can load. Replace "every canonical
> even-k output" with that condition. Byte-order selection is appropriate for
> the validated ASCII alphabet, as specified by the [R radix-sort
> documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/sort.html).
> Keep fixed expected cases such as `AA` versus `TT` and `AACG` versus `CGTT`
> in ordinary tests too. A Danish-locale integration case complements those
> expectations but may be skipped on CI. Restore `LC_COLLATE` even when the
> locale is unavailable or an assertion fails. Host locale does not establish
> which locale a particular user's R session uses.

Evidence:

- `R/ref_kmers.R:515-526` picks the canonical even-k label with
  `min(motif, reverse_complement)`. R compares strings with the current
  collation locale, not by bytes.
- Danish collation treats the pair `"AA"` as a single letter (A-ring), which
  sorts after `"Z"` and therefore after `C`, `G`, and `T`. Labels never
  contain `Z`, but that does not matter. Any comparison whose first
  difference lands on an `"AA"` goes the wrong way. Probe:
  `min("AACG", "CGTT")` is `"AACG"` under `C` and `"CGTT"` under
  `da_DK.UTF-8`, because `"AACG"` compares as A-ring followed by `"CG"`.
- End to end: a canonical store with `AACG`, `ACGT`, and `CCGG` loads under
  `C` and fails under `da_DK.UTF-8` with
  `canonical reference k-mer motif label 'AACG' at motif_idx 1 should be 'CGTT'`.
- Your macOS locale is `da_DK`. R sessions that use the system locale get
  Danish collation (*inference* for RStudio). The probes in this review ran
  with `LANG=C.UTF-8`.
- testthat edition 3 forces `C` collation during tests
  (`testthat::local_reproducible_output()` calls
  `withr::local_collate("C")`), so the test suite cannot catch this.

Impact:

Canonical labels that R computes differently from Rust under `da_DK.UTF-8`,
counted over full k-mer sets:

| k | Canonical labels | Wrong under Danish collation | Example |
| --- | --- | --- | --- |
| 2 | 10 | 1 | `AA` (versus `TT`) |
| 4 | 136 | 18 | `AACG` (versus `CGTT`) |
| 6 | 2,080 | 298 | |
| 8 | 32,896 | 4,806 | |

About 15% of labels are wrong, and one wrong label is enough to fail the load.
In practice every canonical even-k output fails to load for users with Danish
collation. Norwegian collation has the same `"AA"` rule (*inference*, not
probed). Odd k uses the middle-base rule and is not affected.

This is the only place in the R package where locale ordering changes a
result. The other `sort()` calls only build error messages (see Checked And
Not Findings). Python compares strings by code point and is not affected.

Fix:

Compare in byte order with radix ordering, which sorts strings in byte order
whatever the locale is. In `cf_canonical_ref_kmer()`
(`R/ref_kmers.R:517-519`):

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

Why not set the locale before sorting instead:

- `Sys.setlocale()` changes the collation for the user's whole R session, not
  just the package.
- It must be restored with `on.exit()` so the user's setting survives errors.
- `withr::local_collate()` would handle that but adds a dependency.
- Radix ordering is one argument and only affects the call it is used in.

Test:

testthat resets collation to `C`, so the test must set `LC_COLLATE` to
`da_DK.UTF-8` inside the test itself and restore it afterwards. Skip the test
when that locale is not installed. A canonical fixture containing `AACG` is
enough to fail the current code.

### R-02 - High - `end_motif_data_frame(densify = TRUE)` densifies the whole store

> **Review note - two distinct fixes:** Calling the selected-index helper
> below fixes whole-store densification of sparse input. Its dense branch
> still calls `cf_read_array()` before subsetting (`R/end_motifs.R:632-638`),
> so it does not fix whole-store reads for dense input. Ref-kmers also read
> their entire dense array in
> `cf_ref_kmer_dense_frequency_matrix_for_indices()` (`R/ref_kmers.R:991-997`).
> Passing selectors into that helper does not make the underlying read selective.
> Cover dense counts and frequencies as a separate read-path fix.
>
> Do not pass arbitrary index vectors directly to the current CRAN reader and
> assume orthogonal selection. In installed `zarr` 0.4.0,
> `zarr_array$read(selection=...)` takes `min` and `max` on each axis and reads
> that bounding range. It ignores vector order and gaps. Read bounded ranges
> and select relative indices, or use suitable repeated bounded reads, while
> preserving caller order and length-one dimensions. This also constrains R-10.
> Correction still deliberately needs all sample motifs for corrected
> frequencies, as explained in P-08. A motif selector alone does not bound
> its intermediate motif axis.

Evidence:

- `R/end_motifs.R:1341-1345` calls `dense_counts_matrix(x, allow_densify =
  densify)` without the selected rows and motifs, then subsets the result.
- Reference correction goes through the same function. That covers corrected
  data frames with `densify = TRUE` and every `dense_corrected_counts_matrix()`
  call (`R/reference_correction.R:277-283`, `436`).
- The ref-kmer version passes the selection
  (`R/ref_kmers.R:1333-1338`), so this is a bug and not a design choice.

Measured on a sparse store with 50k windows x 1024 motifs (1M stored entries):

- `end_motif_data_frame(window_idxs = 1, motifs = <2 motifs>, densify = TRUE)`:
  about 405 MB extra peak memory. The full dense matrix is 391 MB.
- The same selection through `dense_counts_matrix(..., window_idxs =,
  motifs =)`: 45 MB.

Impact:

The README example at `README.md:112-118` does exactly this. With genome-wide
windows and larger motif sets, a request for a handful of rows runs out of
memory. *Inference*: 290k windows x 1024 motifs is 2.4 GB. A 4+4 motif axis
is far beyond any laptop.

For dense stores, the same line reads the whole counts array on every call
(`R/end_motifs.R:633-637`, and `533` for `sparse_counts_matrix()` on dense
stores). That is smaller but has the same cause.

Fix:

Call `cf_dense_end_motif_matrix_for_indices(x, row_indices, motif_indices,
densify)` instead.

### R-03 - High - `read_lengths()` fails when every group name is numeric

Same bug as Python P-01.

Evidence:

- `R/lengths.R:454` and `457` call `data.table::fread()` without column
  types, so it guesses.
- `R/lengths.R:459-461` converts only `chrom` back to text.
- `R/lengths.R:595` rejects a `group_name` column that is not character.

Probe: group names `1`, `2` and `001`, `002` both fail with
`group_name must contain non-empty character strings`.

Fix:

Read `chrom` and `group_name` as character with `colClasses`. Pass only the
columns that exist in the file, because `fread()` warns about absent names.
Add a test with all-numeric group names.

### R-04 - Medium - Reading `.tsv.zst` needs a separately installed `zstd` program

> **Review note - acceptance criteria:** The `qs2` replacement must support
> streaming frames without a declared decompressed size, including actual
> Rust-writer output, and report corrupted or truncated input clearly. Ensure
> temporary files are removed on success and error. Declare a minimum `qs2`
> version if the file API is not present in every supported release.
> Replace the executable-dependent fixture construction too, so compressed
> loading is exercised when no `zstd` program is installed. The current test
> conditionally tests a valid compressed file and also checks the missing-tool
> error. Its entire test case does not require the executable.

Evidence:

- `R/lengths.R:446-454` calls the `zstd` command-line program and stops if it
  is missing. `.tsv.zst` is the default `cfdna lengths` output.
- The requirement is not mentioned in `README.md`, `DESCRIPTION` (no
  `SystemRequirements`), or `website/docs/guides/fragment_lengths_guide.md`.
- `qs2` is already a required dependency and exports
  `zstd_decompress_file()`. Probe: it decodes a zstd file written in streaming
  mode (no content size in the header, like the Rust writer's
  `zstd::Encoder` at `src/shared/writers.rs:21`). `qs2::zstd_decompress_raw()`
  fails on the same file with `negative length vectors are not allowed`, so
  use the file version.

Impact:

*Inference*: macOS and Windows do not ship `zstd`, so many R users cannot read
lengths output without installing an extra program first.

Fix:

Decompress with `qs2::zstd_decompress_file()` to a temporary file, then
`fread()` it. Drop the command-line path.

### R-05 - Medium - README and website examples that fail or name the wrong flag

> **Review note - cross-language and example scope:** P-02's claim that the
> website ran `ends --all-motifs` is incorrect. It runs `ref-kmers
> --all-motifs`, so both website matrix examples need correction for default
> sparse sample output. Prefer sparse corrected matrices in the general
> recipes. The mode-dependent lengths example issue also affects the Python
> block at `fragment_lengths_guide.md:119-145`, which mixes group and window
> selectors on the same object. Label these as alternatives for their
> respective file modes, rather than suggesting every call can run on the
> same file.

- `README.md:235-239` and `website/docs/guides/fragment_ends_guide.md:301-305`
  call `dense_corrected_counts_matrix(...)` without `allow_densify = TRUE`.
  This errors on default sparse `cfdna ends` output
  (`R/reference_correction.R:420-426`). Same as Python P-02.
- `README.md:80` says `"size"` rows come from `--window-size`. The flag is
  `--by-size` (`src/commands/cli_common.rs:408`).
- `website/docs/guides/fragment_lengths_guide.md:158-163` calls
  `window_metadata(lengths)` and `group_metadata(lengths)` one after the
  other. Only one of them exists for a given file, so running the block as
  written stops with an error.
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

Measured with 256 motifs and full support:

| Windows | Mode | R | Python |
| --- | --- | --- | --- |
| 3k | joint | 6.7 s | 2.2 s |
| 3k | split | 12.4 s | not measured |
| 3k | outside | 10.4 s | not measured |
| 20k | joint | 60 s, 3.5 GB peak process memory | 20 s, about 1.5 GB extra |

*Inference*: genome-wide 10 kb windows (about 290k rows) would take about
15 minutes for `joint` and need far more memory than a laptop has.

Fix:

Same approach as Python P-08, using `Matrix`:

- Line the reference frequencies up with the sample's motif columns once.
- Count supported motifs per row (`frequency > 0`).
- Divide the counts (or the stored sparse values) by
  `frequency * supported_motif_count`.
- Get side frequencies by multiplying the reference matrix with a 0/1 matrix
  that maps each k-mer to its side label.

Build the data frames from the result.

> **Review note - share the acceptance contract with P-08:** Its annotation
> lists the normalization, sparse support, unsupported-count, row matching,
> and frequency invariants that this rewrite must preserve. Use those for
> both packages, with R's one-based public indices. Parse motif sides once
> per unique label. Avoid a dense row-by-full-k-mer intermediate for sparse
> references. Separate the bounded R-02 fix from this larger correction
> refactor so their memory effects can be assessed independently. The table
> compares total R process peak with additional Python memory, so it does not
> establish the ratio of their memory use.

## API And Documentation

### R-07 - Medium - Count matrices have no row or column names

> **Review note - names do not align matrices automatically:** Ordinary matrix
> `rbind()` and `cbind()` bind by position. Adding dimnames helps inspection
> and explicit reindexing but does not prevent silent sample misalignment.
> See the distinction between matrices and data frames in the [R binding
> documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/cbind.html).
> Document explicit matching to a shared motif order before binding.
> Coordinates alone also need not uniquely identify BED rows. Preserve the
> public window/group index in selected-row metadata, especially after
> blacklist filtering, instead of relying solely on `chrom:start-end` names.
> Name raw and corrected columns in every mode, including empty results.
> Existing `unname()` value assertions can remain with separate assertions
> for exact dimnames and order. Removing them alone would not add coverage.

Evidence:

- `sparse_counts_matrix()` and `dense_counts_matrix()` for end motifs and
  ref-kmers return matrices without dimnames (`R/end_motifs.R:530-549`,
  `626-646`, `R/ref_kmers.R:780-797`, `858-881`). Probe on the downstream
  sparse windowed fixture confirms `NULL` dimnames.
- Corrected matrices get column names only in `outside` and `inside` modes
  (`R/reference_correction.R:446-448`, `547-551`).
- Length count matrices have `count_*` column names (`R/lengths.R:639`).
- `dense_counts_vector()` and `dense_frequencies_vector()` are named by motif
  (`R/end_motifs.R:975`, `R/ref_kmers.R:1094`, `1106`).

Impact:

After `motifs = c("AC_GT", ...)` the user gets a matrix that does not say which
column is which motif. In R, dimnames are where users expect labels, and
`cbind()`, `rbind()`, and heatmap functions use them. With motif axes that
differ between samples (P-03), unlabeled columns make silent misalignment
easier.

Fix:

Always set column names to the motif labels. Set row names to window labels
(for example `chrom:start-end`) or group names. `midpoint_array()` and
`profile_array()` could carry dimnames too (group names, length bins, position
starts). The tests currently strip names with `unname()` before comparing
corrected matrices (for example
`downstream_tests/R/test_cfdnalab_r_package.R:198`, `217`, `298`, `327`), so
they would need updating.

### R-08 - Low - Smoothed midpoint values of exactly -1 become NA in R

> **Review note - raise priority and revise the proposed fix:** The installed
> CRAN `zarr` 0.4.0 code masks floating values using `.near(data, fill)`, not
> equality. Its tolerance is `max(eps * max(abs(data), abs(fill)), 1e-12)`,
> with `eps = 1.490116e-08`. Thus loss is not restricted to exact `-1.0` and
> depends on the largest value in the returned slice. This follows from
> `get("zarr_array", asNamespace("zarr"))$public_methods$read` and
> `get(".near", asNamespace("zarr"))`.
> Exact `-1` is reachable with ordinary integer midpoint counts. For
> `savgol=5`, a neighborhood `[23, 0, 2, 0, 0]` gives
> `(-3*23 + 17*2) / 35 = -1` using the coefficients in
> `src/commands/midpoints/smoothing.rs:107-116`.
>
> **Do not use `-Inf`.** For finite data and an infinite fill, both sides of
> the comparison become `Inf`, so every finite value satisfies it and becomes
> `NA`. Use a NaN fill for newly written midpoint arrays. The reader has an
> explicit `is.nan(fill)` branch that masks only NaNs. Update the schema spec,
> writer comment, and writer-backed downstream coverage together. Read actual
> metadata and retain finite negative values, including `-1`, under the new
> fill. Test full and sliced reads with a shortened boundary chunk.
>
> Also audit the other floating arrays. By the same formula, a slice containing
> `2^27` has tolerance 2, so valid zero counts can match a `-1` fill even for
> nonnegative end counts. This is a source-derived consequence, not a new
> measurement. Existing R reads with a finite sentinel cannot reconstruct
> values already masked to `NA`. Treat this as correctness work before
> validation, not a low-priority documentation fix.

Evidence:

- The CRAN `zarr` reader returns `NA` for any stored value equal to the
  array's fill value. Probe: float32 values `[0.5, -1.0, -0.25, 2.0]` with fill
  value `-1` read back as `[0.5, NA, -0.25, 2.0]`. `read()` has no option to
  turn this off.
- The Rust writer chooses fill values outside each array's data domain for
  this reason (`.AI/docs/specs/zarr_schema.md:21-29`).
- `src/shared/zarr.rs:49-51` describes the midpoint fill value as being for
  "non-negative float32 count arrays". Midpoint counts are smoothed by default
  and can be negative (P-04), so `-1.0` is inside their data domain.

Impact:

Rare, because it needs an exact f32 `-1.0`. But it is silent and R-only.
Python returns the real value.

Fix (Rust side):

Use `NaN` or `-Inf` as the fill value for midpoint counts and update the
comment.

### R-09 - Low - Committed man pages are out of date

Evidence:

Running `roxygen2::roxygenise()` (8.0.0, the version in `DESCRIPTION`) in a
scratch copy changes three files:

- `man/dense_corrected_counts_matrix.Rd`
- `man/end_motif_data_frame.Rd`
- `man/sparse_corrected_counts_matrix.Rd`

Each gains the paragraph saying that correction requires
`--orientation both`, is approximate, and needs windows of at least a few
kilobases. `website/docs/generated/loaders/r-api.mdx` is generated from the Rd
files and does not mention `--orientation both`. `NAMESPACE` is in sync.

Fix:

Regenerate the man pages and the website pages. Consider a check that fails
when `roxygenise()` changes tracked files.

### R-10 - Medium - No way to get a selected 3D midpoint array

> **Review note - preserve selectors rather than bounding ranges:** Follow
> R-02's backend-selection constraint. A single read for groups `c(5L, 2L)`
> would read groups 2 through 5 in file order, rather than those two groups
> in requested order. Cover noncontiguous/reversed group and length-bin
> selections, scalar axes, and empty selections without making an invalid
> backend read. Use the existing asymmetric count fixtures to assert values
> as well as shape. `profile_array()` is explicitly a scalar-profile API, so
> its singular selector names are not a broken workflow by themselves.

Evidence:

- `midpoint_array()` reads the whole store (`R/midpoints.R:395-398`).
- `profile_array()` returns one profile (`R/midpoints.R:190-210`).
- Python's `counts_array()` takes the same group and length-bin selectors as
  its data frame method.

Impact:

A heatmap of 50 selected groups means looping over `profile_array()` or
reading everything. For a large grouped store, reading everything is not an
option (*inference* on sizes).

Fix:

Give `midpoint_array()` the selectors from `midpoint_data_frame()`. One read
for a selection is fast: 10 groups x 20 length bins x 2001 positions took
0.01 s.

Related, Low: `profile_array()` takes singular `group` and `length`, while
`midpoint_data_frame()` takes `groups` and `with_lengths`.
`profile_array(m, groups = "LYL1")` fails with `Unused argument(s): groups`.

### R-11 - Shared with the Python review

These apply to R as well. Details are in `python_loader_review.md`.

- P-03: Motif axes differ between samples. `motifs()` and `motif_idx` are
  local to each file.
- P-04: `read_midpoints()` (`R/midpoints.R:19-35`) does not explain values,
  smoothing, or positions either.
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
  open zero-length arrays (`R/end_motifs.R:226-235`).
- P-26: Blacklist rounding in lengths TSVs.
- P-27: Complement collapse is not reliably recorded, and correction can use
  the wrong reference denominator. Its annotation applies to R too.
- P-29: End-motif sparse validation lacks the sorted/unique coordinate check
  present in the ref-kmer reader. This can be fixed without a shared-module
  refactor.

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

## Tests

- No test runs under a non-`C` collation, so R-01 cannot be caught. testthat
  forces `C` collation, so the test must set the locale itself.
- No test uses all-numeric group names (R-03).
- No test would catch R-02, because the result is correct. A sparse fixture
  with a large declared shape (for example 100000 x 1000) and a few stored
  entries would make a full densification allocate about 800 MB, so a
  regression would show up as a failure.
- The `.tsv.zst` test depends on the `zstd` program being installed
  (`tests/testthat/test-lengths.R:511-527`).
- Corrected matrix comparisons strip names with `unname()`, so the dimnames
  gaps in R-07 are built into the expectations.
- Motif-group correction values are untested (see R-11).

> **Review note - avoid an allocation-failure regression test:** A declared
> `100000 x 1000` shape does not guarantee failure. Allocating 800 MB may
> succeed, or may kill the entire R process. Keep ordinary fixtures tiny and
> assert selected values, order, and retained empty rows. Add a small fake
> read/conversion boundary that rejects whole-store materialization for a
> bounded selection, rather than depending on allocator behavior or elapsed
> time. This should distinguish the sparse densification fix from the dense
> selective-read fix in R-02.
