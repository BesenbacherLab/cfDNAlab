# Output Loader Obvious Fixes Plan

Status: proposed first implementation batch. No fixes are implemented by this document.

Implement small, compatible fixes from the [Python review](../major_reviews/python_loader_review.md)
and [R review](../major_reviews/r_loader_review.md). Keep larger performance work
and unresolved API or scientific choices separate.

## Scope And Agreed Constraints

- Continue reading previously supported, valid cfDNAlab outputs. Do not require
  regenerated outputs, new mandatory metadata, or a schema migration in this batch.
- Defer any change that requires breaking compatibility for a separate decision.
- Keep existing public names, data frame columns, matrix shapes, index bases,
  and selector ordering.
- Do not add a shared motif-axis API or raw end-motif frequencies.
- Leave reference-correction support, unsupported-count handling, and frequency
  normalization unchanged. The earlier selected-only validation proposal is
  deferred pending discussion, including for count-only outputs.
- CLI format decisions govern the readers. Do not change floating fill values
  to accommodate the R backend in this batch.
- Do not introduce a generic loader framework or rewrite the correction engine.

## Invariants

Before changes to masking, aggregation, or normalization, record the relevant
invariants from this list beside the implementation work:

- Labels retain their exact text and identify the same rows as before.
- Python selectors stay zero-based. Public R selectors and index columns stay
  one-based. Genomic and profile coordinates keep their existing meanings.
- Selection preserves requested order and length-one dimensions.
- Sparse and dense extraction return the same values for the same selection.
- Finite negative midpoint values remain valid. Counts from ends and lengths
  remain finite and nonnegative.
- Reference correction formulas, full-reference support, side aggregation,
  corrected-frequency denominators, and unsupported-count policies stay unchanged
  in this batch.
- Sparse reads and validation must not introduce whole-store densification.

## 1. Preserve Length-Table Labels

Findings: P-01, R-03.

Files: `py-cfdnalab/src/cfdnalab/lengths.py`, `r-cfdnalab/R/lengths.R`.

- Parse `chrom` and `group_name` as text before numeric inference changes them.
- Remove the casts that attempt to recover chromosome text after parsing.
- Preserve literal labels such as `NA`. In R, check the parser's missing-value
  rules as well as `colClasses`. Keep numeric-column validation intact.
- Apply the same parsing rules to plain and Zstandard-compressed input.

Acceptance cases: `001` beside `1`, `1.10` beside `1.1`, all-numeric group names,
literal `NA`, and nonnumeric names. Assert exact metadata and group-name
selection results. Use legal fragment length columns such as `count_30`.

## 2. Fix R Canonical Label Comparison

Finding: R-01.

File: `r-cfdnalab/R/ref_kmers.R`.

- Use explicit byte-order comparison for even-k canonical labels, matching Rust.
- Keep the odd-k middle-base rule.
- Do not change the user's session collation in production code.
- If tests use `withr` directly, declare it in `Suggests`.

Acceptance cases: fixed expected representatives for `AA` versus `TT` and
`AACG` versus `CGTT`, plus a canonical store containing an affected label.
The Danish-locale case must confirm that the requested collation actually took
effect, restore it on exit, and skip when unavailable. Keep ordinary fixed-label
cases too so not all coverage depends on an installed locale.

## 3. Make Contained Selection And Error-Handling Fixes

Findings: the sparse part of R-02, P-06, parts of P-09, P-20, P-28.

Files: Python `ends.py`, `ref_kmers.py`, and `midpoints.py`, plus the matching
R end-motif reader.

- In R, pass selected row and motif indices into the existing helper before
  densifying sparse end-motif counts.
- In Python, replace individual midpoint-profile reads with a bulk orthogonal
  read. Keep all three dimensions and requested order.
- For Python dense-to-sparse extraction, select the dense rows and columns
  before conversion instead of reading the full dense matrix.
- Vectorize the Python ref-kmer sparse-coordinate order check.
- Add a separate Python midpoint test with writer-consistent coordinates and
  distinguishable values per axis. Leave the existing negative-position test
  and its fixture unchanged.
- Narrow missing-array exception handling so corrupt metadata and I/O failures
  retain useful context. Do not report every exception as a missing array.

Acceptance cases: scalar, reversed, noncontiguous, and empty selections.
Verify values and selected-axis identity, and distinguish missing arrays from
metadata or I/O errors.

Retain existing validation without adding new scans. R already checks sparse
index integer values, bounds, and finite nonnegative counts. Its current sparse
selection uses membership and matching rather than assuming coordinate order.
Adding end-motif order and uniqueness checks or broader malformed-store
validation is outside this batch. Vectorizing an existing Python check retains
its behavior while removing the per-entry Python loop.

Keep these limits explicit:

- R's selected-index helper still reads whole dense arrays. The sparse fix does
  not resolve the dense-read part of R-02.
- R backend selections read bounding ranges. Selecting distant indices can
  therefore still materialize the whole store. A complete dense-read fix needs
  separate work with bounded reads.
- Vectorized masks still scan all sparse entries. Repeated-access indexing is
  deferred, and P-07 must remain open for that part.

Use small value regressions and a fake read boundary where useful. Do not rely
on allocation failure, timing limits, or forbidding a particular helper name.

## 4. Simplify Compressed-Table Reading And Preserve Errors

Findings: R-04 and the compressed-input part of P-28.

The current R path works when the `zstd` command-line program is installed.
This is a portability improvement, not a claim that it fails in that setup.

- Evaluate `qs2::zstd_decompress_file(input_path, temporary_path)`, followed by
  `data.table::fread()`. `qs2` is already an R dependency.
- Require support for actual Rust-written streaming frames, clear errors for
  corrupt or truncated input, and cleanup on success and error.
- Pin a supported minimum version after checking those requirements. The file
  API was introduced in 0.1.6 according to the
  [upstream changelog](https://github.com/qsbase/qs2/blob/main/ChangeLog).
- Use the file API for new compressed test fixtures if the replacement is
  adopted. Leave existing tests and fixture construction unchanged unless the
  user approves changing them.
- Remove Python's redundant command-line fallback and include the original
  parsing or decompression error in the user-facing message.

Claude's review reports successful streaming-input checks with `qs2` 0.1.7.
This plan does not treat that as independent verification of the replacement.
If the file API does not meet the requirements, keep the working path and defer
the replacement rather than silently weakening error handling.

## 5. Correct Documentation And Blacklist Precision

Findings: P-02 to P-05 documentation, P-16, P-18, P-19, P-26, R-05, R-09,
and the unused `LazyData` setting in R-13.

- Use sparse corrected matrices in general examples. Keep any dense example
  explicitly opted in and limited to a small result.
- Explain stored counts, optional weighting, fractional overlap assignment,
  midpoint smoothing, per-bp bin averages, interval-relative orientation, and
  normalization only when `eligible_intervals > 0`.
- Explain motif direction and storage. Document that motif indices are local
  to each file and samples need explicit alignment by labels.
- Separate global, windowed, and grouped examples so each block is valid for
  its stated mode.
- Complete missing Python argument and return documentation, including shapes,
  selector order, blacklist filtering, and densification behavior.
- Preserve round-trip precision for `blacklisted_fraction` in new lengths TSVs.
  Keep reading existing rounded files and leave count rounding unchanged.
  Document that discarded digits in old files cannot be recovered.
- Remove the unused R `LazyData` setting.
- Update roxygen sources and regenerate man pages. Then regenerate website
  loader documentation. Do not edit generated files by hand.

Keep broader documentation restructuring and new analysis recipes out of this
batch. Update current specs only when the corresponding behavior is implemented.

## Verification And Completion

- Add new test cases for fixes. Do not modify or remove existing tests, their
  assertions, or their fixtures without the user's explicit permission.
- Existing tests may fail after a change. Report and investigate the failure
  rather than changing the test to make it pass. Request permission if the
  intended behavior requires updating an existing test.
- Write small regressions with hand-derived expectations before each fix.
  Re-read the relevant code and existing tests at implementation time.
- Cover both current and previously supported valid output fixtures.
- Add meaningful assertions for values and selected-axis identity. Do not
  duplicate the implementation to generate expected answers.
- The agent does not run tests. The user runs Python, R, Rust, and downstream
  suites. Provide the relevant commands when handing over implementation work.
- After Rust code changes, run `cargo check`. After Rust test changes, run
  `cargo check --tests --features testing`. Include non-default features touched
  by a patch. Run `cargo check --all-features` after a major refactor.
- Perform appropriate Python/R syntax and documentation checks. Do not treat
  compilation or static checks as proof that regressions pass.
- Mark only completed portions of review findings resolved. Record remaining
  limitations, especially R-08, dense R reads, and repeated sparse scans.

## Deferred Engineering Work

- P-08/R-06: replace long-table reference correction with bounded matrix or
  sparse operations after the correction-selection policy is resolved. Do not
  include the earlier selected-only proposal as an agreed requirement.
- R-08: investigate the R backend's floating-value masking separately. No CLI
  fill-value change is proposed. A reader adaptation must preserve encoded
  values and absent-chunk handling without changing global tolerances or
  guessing replacements for missing values. The mutable fill-value binding is
  only a candidate identified from source inspection, not a verified fix.
- P-29: defer additional malformed-store validation. Add checks only where a
  concrete failure or an algorithm's requirements justify their load-time and
  memory cost. Do not add duplicate R checks or validate again on every access.
- Complete selective dense reading in R and repeated sparse row indexing in
  both packages. Avoid a large cache or backend abstraction in this first batch.
- Further large-axis loading optimizations and documentation deduplication.

## Focused Decisions For Later Discussion

| Topic | Decision needed |
| --- | --- |
| Shared motif axes | Cohort/model use cases, fixed panels, zero versus uncounted motifs, and necessary target metadata |
| Raw end-motif frequencies | Which counted subspace defines the denominator |
| Reference-correction selection | Whether count-only outputs may ignore unsupported positive counts outside the selection. Discuss correction factors and corrected-frequency denominators separately, including side aggregation |
| API consistency | Python index-column renames, R matrix names, and how filtered arrays expose row identity |
| Midpoint conveniences | Center-relative coordinates, extra bin metadata, selected R arrays, and settings accessors |
| Empty sparse outputs | Whether metadata-only loading is useful and what selection/correction should mean |
| Complement collapse | Whether to retain the option and how correction should handle older outputs with unknown collapse status |

P-27 remains deferred because adding a flag to future stores does not resolve
the missing information in existing stores. Its one-base example is also
misleading under `orientation="both"`, where normalization can cancel out.
Any later regression should use longer motifs whose complement frequencies
can differ despite reverse-complement symmetry.

Required new metadata, schema changes, and compatibility breaks need a separate
decision before implementation.

## Context For Deferred Items

The `-1` midpoint example establishes possibility, not prevalence. For the
allowed `savgol=5` setting, coefficients `[-3, 12, 17, 12, -3] / 35` applied to
integer counts `[23, 0, 2, 0, 0]` give `-1` at the center. This is a constructed
case, not evidence that exact `-1` occurs often with the default 165 bp window.
The installed R backend also masks approximately equal values with a tolerance
that grows with the largest magnitude in the read slice. Exact equality is
therefore not the only issue. Neither observation determines the CLI fill-value
policy. Valid frequency arrays are bounded by `[0, 1]`, so large-count masking
examples do not establish the same failure for valid frequencies.

The current Python and R count-only correction methods call the data frame
correction path, which also calculates frequencies. That path corrects the full
sample motif axis, calculates frequencies, and only then selects output motifs.
A count-only result and a frequency result have different dependencies, but
separating their implementations must not silently change the scientific or
unsupported-count rules. Keep current behavior until that policy is discussed.
