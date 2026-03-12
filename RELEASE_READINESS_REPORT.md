# cfDNAlab Release Readiness Report

**Date:** 2026-03-12
**Version reviewed:** 0.1.0
**Branch:** `main`
**Reviewed commands:** All 9 commands listed in the root README (`fcoverage`, `midpoints`, `lengths`, `gc-bias`, `ref-gc-bias`, `coverage-weights`, `bam-to-bam`, `bam-to-frag`, `frag-to-bam`) and their imported shared utilities.

---

## Executive Summary

| Area | Status |
|------|--------|
| **Build** | Clean — compiles with 0 errors (debug and release) |
| **Tests** | All 349 tests pass, 0 failures, 4 ignored doctests |
| **Clippy** | **425 warnings** (214 auto-fixable) |
| **TODOs in source** | **~65 TODO/FIXME** comments in released command paths |
| **License** | **Missing** — no LICENSE file |
| **Cargo.toml metadata** | **Missing** license, description, repository, authors, keywords |
| **Changelog** | Version mismatch (`0.0.1` vs Cargo.toml `0.1.0`) |

---

## TODO Checklist

### TIER 1 — Must Fix Before Release

These are blockers. The tool should not be publicly released until these are resolved.

---

- [ ] **T1-1: Add a LICENSE file**
  - **Location:** Repository root (missing)
  - **Severity:** BLOCKER
  - **Description:** There is no `LICENSE` or `LICENCE` file anywhere in the repository. Without a license, users have no legal clarity on usage rights. Most package registries and institutions require a license for adoption. Choose an appropriate open-source license (e.g., MIT, Apache-2.0, GPL) and add the file to the repo root.

---

- [ ] **T1-2: Add required Cargo.toml metadata**
  - **Location:** `Cargo.toml` (lines 1–4)
  - **Severity:** BLOCKER
  - **Description:** The `[package]` section is missing all standard metadata fields expected for a published crate:
    - `license` (or `license-file`)
    - `description`
    - `repository`
    - `authors`
    - `readme`
    - `keywords`
    - `categories`
    - `homepage` (optional)
    - `documentation` (optional)
  - These are required for crates.io publishing and expected by users evaluating the tool.

---

- [ ] **T1-3: Validate or rewrite `gc_bias/interpolation.rs` (auto-generated code)**
  - **Location:** `src/commands/gc_bias/interpolation.rs:1`
  - **Severity:** BLOCKER
  - **Description:** Line 1 reads: `///! NOTE: This code was generated. TODO: Validate that it's correct.` This file implements polynomial interpolation for GC bias correction — a core feature of the tool. The entire module (~500 lines) was LLM-generated and explicitly marked as unvalidated. Before release, the code must be:
    1. Manually reviewed for correctness (especially the Gauss-Jordan solver at `solve_sym_posdef`, the monotonic enforcement at `enforce_monotonic_segment`, and the weighted polynomial fitting)
    2. Validated against known-good reference implementations or analytical test cases
    3. The warning comment removed once validated
  - The existing tests in `tests/test_interpolation.rs` cover basic scenarios but may not be sufficient for full validation.

---

- [ ] **T1-4: Validate or rewrite `shared/frag_file.rs` (auto-generated code)**
  - **Location:** `src/shared/frag_file.rs:1`
  - **Severity:** BLOCKER
  - **Description:** Line 1 reads: `// TODO: Fix this. Just generated but chatty doesn't know frag files (finaledb) so it invents stuff`. This module is used by `frag-to-bam` to parse fragment files. The comment explicitly states the LLM "invents stuff" about the format. The code must be validated against the actual finaledb/frag file specification and the comment removed.

---

- [ ] **T1-5: Clean up README placeholder TODOs and missing references**
  - **Location:** `README.md` (multiple lines)
  - **Severity:** BLOCKER
  - **Description:** The public-facing README contains unfinished placeholders that make the project look incomplete:
    - **Line 9:** States `"The package is in alpha-stage (being developed)"` — update for release
    - **Line 73:** `[TODO: Not that simple]` and `(TODO on samtools!)` in the FAQ
    - **Line 280:** `[REFS]` — missing citation for fragment length cancer detection studies
    - **Line 313:** `[REFS]` — missing citation for midpoint coverage studies
    - **Line 340:** `[TODO: Note on how to get griffin-like profiles]`
    - **Line 348:** `[TODO: Add output-prefix for remaining commands]`
    - **Line 437:** `[TODO: Check correctness]` for file column documentation
  - All of these should be resolved (filled in or removed) before release.

---

- [ ] **T1-6: Fix changelog version mismatch**
  - **Location:** `CHANGELOG` (line 3) vs `Cargo.toml` (line 3)
  - **Severity:** HIGH
  - **Description:** `CHANGELOG` says `v/ 0.0.1` while `Cargo.toml` says `0.1.0`. These must be consistent. Additionally, the `v/` formatting looks like a typo (should be `v0.1.0` or similar).

---

### TIER 2 — Strongly Recommended Before Release

These are not strict blockers but represent significant quality issues. Shipping without addressing them creates risk.

---

- [ ] **T2-1: Fix typo in user-facing error message (`indel_mode.rs`)**
  - **Location:** `src/shared/indel_mode.rs:30`
  - **Severity:** HIGH (user-facing)
  - **Description:** The error message reads `"Use 'ignore', 'adjus', or 'skip'"` — missing the `t` in `'adjust'`. This is shown to users when they provide an invalid `--indel-mode` value. Simple one-character fix.

---

- [ ] **T2-2: Address 425 clippy warnings**
  - **Location:** Entire codebase
  - **Severity:** HIGH
  - **Description:** `cargo clippy` reports 425 warnings. Breakdown:
    | Warning type | Count | Auto-fixable |
    |-------------|-------|-------------|
    | Doc list item overindented | 131 | Yes |
    | Unnecessary reference creation | 25 | Yes |
    | Loop variable type hints | 15 | Yes |
    | Borrowed expression implements required traits | 15 | Yes |
    | Casting to the same type (u32→u32, usize→usize, etc.) | 37 | Yes |
    | Module has same name as containing module | 13 | No (rename) |
    | Very complex type (needs type alias) | 13 | No (refactor) |
    | Redundant field names in struct init | 11 | Yes |
    | Too many arguments (8/7) | 5 | No (refactor) |
    | Manual reimplementation of std methods | 7 | Yes |
    | Other | ~153 | Mixed |
  - **Action:** Run `cargo clippy --fix --lib -p cfdnalab` to auto-fix ~214 warnings, then manually address the remainder. For a public release, zero clippy warnings is the standard expectation.

---

- [ ] **T2-3: Replace `unreachable!()` with proper error handling in released commands**
  - **Locations:**
    - `src/commands/bam_to_bam/bam_to_bam.rs:352`
    - `src/commands/lengths/lengths.rs:770`
    - `src/commands/midpoints/midpoints.rs:525`
    - `src/commands/fcoverage/fcoverage.rs:426, 433, 479, 513`
    - `src/commands/bam_to_frag/bam_to_frag.rs:389`
  - **Severity:** HIGH
  - **Description:** These are all the same pattern — a match on `(gc_weight_option, correct_gc)` with `(Some(_), false) => unreachable!()`. While logically sound today (if `correct_gc` is false, `gc_weight_option` should be None), the `unreachable!()` macro causes an immediate panic with no useful message if the invariant is ever violated. In a released tool processing large datasets, a panic means lost computation.
  - **Fix:** Replace with `anyhow::bail!("unexpected gc weight when gc correction is disabled")` or restructure to `if correct_gc { ... } else { ... }` which eliminates the unreachable arm entirely.

---

- [ ] **T2-4: Fix NaN panic risk in interpolation sort**
  - **Location:** `src/commands/gc_bias/interpolation.rs:280-281`
  - **Severity:** HIGH
  - **Description:** The code sorts anchor points by distance using:
    ```rust
    left_real.sort_by(|lhs, rhs| lhs.2.partial_cmp(&rhs.2).unwrap());
    right_real.sort_by(|lhs, rhs| lhs.2.partial_cmp(&rhs.2).unwrap());
    ```
    `partial_cmp` returns `None` for NaN values, and `.unwrap()` on `None` panics. If any distance calculation produces NaN (e.g., from degenerate coordinates), the GC bias correction pipeline crashes.
  - **Fix:** Use `.unwrap_or(std::cmp::Ordering::Equal)` or `total_cmp()` (available on f64 in recent Rust editions).

---

- [ ] **T2-5: Add integration tests for `coverage-weights` command**
  - **Location:** `tests/` (missing file)
  - **Severity:** HIGH
  - **Description:** The `coverage-weights` command has **no integration test file at all**. This is a released command that users depend on for genomic smoothing. The striding and triangular overlap logic is tested in isolation, but the full command pipeline (BAM → stride bins → smoothed weights → TSV output) is never tested end-to-end.

---

- [ ] **T2-6: Add integration tests for `midpoints` command**
  - **Location:** `tests/test_profile_groups_command.rs` (partial), `tests/test_heatmap_inputs.rs` (1 test)
  - **Severity:** HIGH
  - **Description:** The midpoints command has no true integration test. `test_profile_groups_command.rs` tests some sub-components, and `test_heatmap_inputs.rs` has exactly 1 test for heatmap rendering. The full pipeline (BAM → grouped interval counting → profile output) is untested. Additionally, `counting_by_group.rs:335` is explicitly marked `// TODO: Test!!` — the `view_ndarray3_group_len_pos()` zero-copy view method has zero test coverage despite being a core output reshaping function.

---

- [ ] **T2-7: Add integration tests for `ref-gc-bias` command**
  - **Location:** `tests/test_ref_gc_bias.rs` (75 lines, 2 tests)
  - **Severity:** MEDIUM-HIGH
  - **Description:** The `ref-gc-bias` command has only 2 unit tests (75 lines total). No integration test exercises the full pipeline (2bit reference → windowed GC counting → output). This command produces reference data that all downstream GC correction depends on — errors here propagate to every sample.

---

- [ ] **T2-8: Validate or re-examine auto-generated test file**
  - **Location:** `tests/test_coverage.rs:3`
  - **Severity:** MEDIUM-HIGH
  - **Description:** Line 3 reads: `// TODO: Check manually - generated but not validated!` This file contains 31 tests for the `Coverage` struct (prefix sums, blacklist masking, window queries). While the tests pass, their correctness has never been manually verified. A test that passes but checks wrong values provides false confidence.

---

- [ ] **T2-9: Validate uncertain GC bias test**
  - **Location:** `tests/test_gc_bias_windows.rs:372`
  - **Severity:** MEDIUM
  - **Description:** Test `merges_crossing_files_and_scales_once_per_window()` contains the comment `// TODO: Validate this`, indicating the test author was uncertain about the expected values. This test exercises cross-tile merging logic which is critical for correctness.

---

- [ ] **T2-10: Remove dead `keep_temp = false` pattern**
  - **Locations:**
    - `src/commands/lengths/lengths.rs:355`
    - `src/commands/midpoints/midpoints.rs:267`
    - `src/commands/fcoverage/fcoverage.rs:523`
    - `src/commands/wps/wps.rs:564`
    - `src/commands/gc_bias/gc_bias.rs:355`
  - **Severity:** MEDIUM
  - **Description:** Five commands contain `let keep_temp = false;` followed by an if/else branch where the `else` (keeping temp files) is dead code. The TODO comment on several says `"Make cli arg behind a feature for dev purposes?"`.
  - **Action:** Either (a) remove the dead branch entirely, or (b) implement it as a `--keep-temp` CLI flag gated behind a dev feature. Dead code in a release is confusing.

---

- [ ] **T2-11: Remove dead `quiet = false` pattern**
  - **Locations:**
    - `src/commands/bam_to_bam/bam_to_bam.rs:51`
    - `src/commands/bam_to_frag/bam_to_frag.rs:54`
  - **Severity:** LOW-MEDIUM
  - **Description:** Same pattern as `keep_temp` — hardcoded `let quiet = false;` with dead `if quiet` branches. Either implement `--quiet` as a CLI flag or remove the dead code.

---

- [ ] **T2-12: Remove commented-out code blocks**
  - **Locations:**
    - `src/shared/fragment_iterator.rs:777-795` — commented-out `fragments_with_records_from_iter` function
    - `src/commands/coverage_weights/striding.rs:150-195` — commented-out `get_overlapping_normalization` function
    - `src/shared/bam.rs:37+` — commented-out `bam_header_contigs_with_len` function
    - `src/commands/ends/ends.rs:9+` — commented-out `run` function
    - `src/lib.rs:7-8` — commented-out `pub use` statements
  - **Severity:** MEDIUM
  - **Description:** Large blocks of commented-out code are confusing to contributors and reviewers. If the code is needed later, it's in git history. Remove for a clean release.

---

- [ ] **T2-13: Remove stale TODO comment on `read.rs`**
  - **Location:** `src/shared/read.rs:75`
  - **Severity:** LOW
  - **Description:** Comment says `// TODO: Requires testing!` for `parse_md_tag()`. However, `tests/test_parse_md_tags.rs` exists and does test this function. The TODO is stale and should be removed.

---

### TIER 3 — Should Fix (Error Handling & Robustness)

These improve the robustness of the tool and prevent panics on edge cases or corrupted input.

---

- [ ] **T3-1: Replace `.expect()` chains in `gc_bias/correct.rs` with error propagation**
  - **Location:** `src/commands/gc_bias/correct.rs:31, 35, 39, 43, 68`
  - **Severity:** HIGH
  - **Description:** `GCCorrector::from_package()` has 4 consecutive `.expect()` calls on edge vector first/last elements, plus one `.expect("fragment end precedes start")`. If a user provides a corrupted or truncated `.npz` GC correction file, these panic instead of returning a helpful error.
  - **Fix:** Replace with `.ok_or_else(|| anyhow!("..."))?` to propagate errors.

---

- [ ] **T3-2: Replace `.expect()` calls in `frag_to_bam.rs` chromosome lookup**
  - **Location:** `src/commands/frag_to_bam/frag_to_bam.rs:191, 226`
  - **Severity:** HIGH
  - **Description:** `expect("chromosome length available for first/next chromosome")` panics if a chromosome from the fragment file doesn't exist in the reference. This can happen if the user provides mismatched files.
  - **Fix:** Return a user-friendly error with the chromosome name.

---

- [ ] **T3-3: Replace `.expect()` in reducer contig lookup**
  - **Location:** `src/commands/fcoverage/fcoverage.rs:494-498`
  - **Severity:** HIGH
  - **Description:** `expect("missing contig length")` panics if a chromosome from temp files doesn't match the contig map. Could happen with corrupted intermediate files.
  - **Fix:** Use `.ok_or_else(|| anyhow!("Chromosome '{}' not found in contig map", chr))?`.

---

- [ ] **T3-4: Replace `unwrap()` calls in reducer K-way merge**
  - **Location:** `src/commands/fcoverage/reducer.rs:263, 274, 281, 507, 517, 526`
  - **Severity:** HIGH
  - **Description:** The core merge loop has 6 unwrap/expect calls that assume heap/stream invariants hold. If any tile produces unexpected output, these panic mid-reduction (after potentially hours of computation).
  - **Fix:** Propagate errors with context about which tile/stream caused the failure.

---

- [ ] **T3-5: Replace `unwrap()` in ndarray slice operations**
  - **Location:** `src/commands/gc_bias/outliers.rs:245-246, 253`
  - **Severity:** MEDIUM-HIGH
  - **Description:** `.as_slice().unwrap()` and `.as_slice_mut().unwrap()` on ndarray arrays. These panic if the array is not contiguous in memory. While typically safe for freshly-allocated arrays, ndarray does not guarantee contiguity after slicing/transposing operations.
  - **Fix:** Check `.is_standard_layout()` first or use `.as_slice().ok_or_else(|| ...)?`.

---

- [ ] **T3-6: Replace mutex `.lock().unwrap()` with poison handling**
  - **Location:** `src/commands/midpoints/counting_by_group.rs:399, 417, 426, 428`
  - **Severity:** MEDIUM-HIGH
  - **Description:** The parallel merge path for midpoint profile groups uses `.lock().unwrap()` and `Arc::try_unwrap().expect()`. If a worker thread panics (e.g., from an I/O error), the mutex gets poisoned and all subsequent `.lock().unwrap()` calls cascade into more panics, making the root cause hard to diagnose.
  - **Fix:** Use `.lock().map_err(|e| anyhow!("worker thread panicked: {}", e))?`.

---

- [ ] **T3-7: Replace `lengths.rs:401` unchecked index with bounds check**
  - **Location:** `src/commands/lengths/lengths.rs:401`
  - **Severity:** MEDIUM-HIGH
  - **Description:** `all_bins[0]` is accessed without checking if `all_bins` is empty. If no windows are processed (e.g., all chromosomes skipped by blacklist), this panics. The same pattern appears at line 415.
  - **Fix:** Add `if all_bins.is_empty() { return Ok(()); }` before the access.

---

- [ ] **T3-8: Guard `heatmap.rs` unwraps on edge vectors**
  - **Location:** `src/shared/plotters/heatmap.rs:358, 359, 546, 547, 964`
  - **Severity:** MEDIUM
  - **Description:** Multiple `.first().unwrap()` and `.last().unwrap()` calls on edge vectors without emptiness checks. If called with empty data, these panic.
  - **Fix:** Add emptiness guards before accessing first/last elements.

---

- [ ] **T3-9: Fix duplicate nested condition in `lengths.rs`**
  - **Location:** `src/commands/lengths/lengths.rs:470-471`
  - **Severity:** LOW-MEDIUM
  - **Description:** Copy-paste artifact:
    ```rust
    if opt.gc.gc_file.is_some() {
        if opt.gc.gc_file.is_some() {  // DUPLICATE CHECK
    ```
    The inner check is redundant. While not a bug, it suggests the surrounding logic may not have been carefully reviewed.

---

- [ ] **T3-10: Add division-by-zero guard in `striding.rs`**
  - **Location:** `src/commands/coverage_weights/striding.rs:146`
  - **Severity:** MEDIUM
  - **Description:** `bins[i].avg_overlap_coverage = sum_cov / (sum_w as f32);` — if `sum_w` is 0, this produces infinity. While the triangular weight loop should always produce at least weight 1, a defensive check prevents silent data corruption.
  - **Fix:** `if sum_w > 0 { sum_cov / (sum_w as f32) } else { 0.0 }`.

---

- [ ] **T3-11: Validate `StrideBin` invariant (start <= end)**
  - **Location:** `src/commands/coverage_weights/striding.rs:24`
  - **Severity:** LOW-MEDIUM
  - **Description:** `size()` uses `self.end.saturating_sub(self.start)` with a TODO comment `"Should not happen?"`. If `start > end` is truly impossible, add a debug assertion at construction time. If it can happen, handle it explicitly rather than silently returning 0.

---

- [ ] **T3-12: Replace `ProgressStyle` template `.unwrap()` calls**
  - **Locations:**
    - `src/commands/lengths/lengths.rs:140`
    - `src/commands/midpoints/midpoints.rs:164`
    - `src/commands/fcoverage/fcoverage.rs:223`
    - `src/commands/ref_gc_bias/ref_gc_bias.rs:136`
  - **Severity:** LOW
  - **Description:** All progress bar templates use `.template("...").unwrap()`. While the hardcoded template strings are valid, any future change could introduce a panic. Consider using `.template("...").expect("hardcoded template")` at minimum for clarity.

---

### TIER 4 — Nice to Have (Code Quality & Maintainability)

These are not urgent but improve the codebase for long-term maintenance.

---

- [ ] **T4-1: Add type aliases for complex iterator return types**
  - **Location:** `src/shared/fragment_iterator.rs` (multiple functions)
  - **Severity:** LOW
  - **Description:** Clippy flags 13 "very complex type" warnings. Functions like `fragments_with_segments_from_bam` return `PairingAdapter<impl Iterator<Item = Result<InputItem<FragmentWithSegments>>>, WithSegmentsPairer, SegmentedReadInfo, FragmentWithSegments>`. Type aliases would make the API more readable.

---

- [ ] **T4-2: Gate always-compiled non-released modules behind features**
  - **Location:** `src/commands/mod.rs:9, 14, 27, 28`
  - **Severity:** LOW
  - **Description:** The `ends`, `fragment_kmers`, `transitions`, and `visualize_positions` modules are compiled unconditionally (no `#[cfg(feature = ...)]` gate), unlike the released commands. This adds to compilation time and binary size even when unused.

---

- [ ] **T4-3: Use `checked_mul` for stride calculations in `counting_by_group.rs`**
  - **Location:** `src/commands/midpoints/counting_by_group.rs:349`
  - **Severity:** LOW
  - **Description:** `let strides = ((p * l) as usize, 1usize, l as usize);` — the multiplication `p * l` could theoretically overflow before the cast. Use `p.checked_mul(l).ok_or_else(|| ...)?` for defensive coding.

---

- [ ] **T4-4: Convert `debug_assert` to `ensure!` in release-critical paths**
  - **Location:** `src/commands/midpoints/midpoints.rs:582-589, 603-610`
  - **Severity:** LOW-MEDIUM
  - **Description:** Window position bounds are validated with `debug_assert!` only — these checks are stripped in release builds. If a window start exceeds `u32::MAX` (cast from u64), the index wraps silently. Consider using `anyhow::ensure!()` for production validation, or at minimum document why the invariant is guaranteed.

---

- [ ] **T4-5: Consider `TryFrom`/`TryInto` for integer casts in hot paths**
  - **Location:** Multiple files (120+ `as usize`/`as u32`/`as u64` casts)
  - **Severity:** LOW
  - **Description:** The codebase has 120+ bare `as` casts between integer types. For human genomes these are always safe (max coordinate ~249M fits in u32), but for non-standard organisms or future use cases, some could silently truncate. High-risk locations:
    - `bam_to_bam.rs:417-418` — `rec.pos() as u32` (i64 → u32)
    - `fcoverage/tiling.rs:336-342` — `abs_start as u32` (u64 → u32)
    - `fcoverage/fcoverage.rs:1099-1100` — blacklist coordinates u64 → u32
  - **Action:** Consider adding a note in the README about supported coordinate ranges, or use `.try_into()` for user-facing inputs.

---

- [ ] **T4-6: Rename cryptic variables in `fcoverage/tiling.rs`**
  - **Location:** `src/commands/fcoverage/tiling.rs:193`
  - **Severity:** LOW
  - **Description:** TODO comment says `"Rename variables to meaningful names"`. Variables like `cs`, `ce`, `k_lo`, `k_hi` make the bin/window calculation hard to follow. Rename to `core_start`, `core_end`, `first_window_idx`, `last_window_idx` or similar.

---

- [ ] **T4-7: Clean up `base.rs` TODO comment**
  - **Location:** `src/shared/base.rs:32`
  - **Severity:** LOW
  - **Description:** Comment says `// TODO: What is "anything else" possibly?` about the nucleotide encoding. "Anything else" maps to code 4 (N) and includes IUPAC ambiguity codes (R, Y, S, W, etc.), soft-masked lowercase, and any non-ACGT byte. Resolve the TODO by documenting this.

---

- [ ] **T4-8: Address `with_records_fragment.rs` clone TODO**
  - **Location:** `src/shared/fragment/with_records_fragment.rs:98`
  - **Severity:** LOW
  - **Description:** TODO says `"Avoid cloning. Would like to keep reusing oriented_pair_from_read_info but perhaps an owned version of it is needed?"`. BAM records are cloned unnecessarily in the `bam-to-bam` pipeline. Not a correctness issue but wastes memory proportional to fragment count.

---

- [ ] **T4-9: Address `blacklist/load.rs` streaming TODO**
  - **Location:** `src/shared/blacklist/load.rs:18-19`
  - **Severity:** LOW
  - **Description:** TODO says `"plumb through a streaming reader so we do not materialise every interval before merging when very large inputs are used."` Currently all blacklist intervals are loaded into memory before merging. For very large blacklist files this could cause high memory usage.

---

- [ ] **T4-10: Document fragment counter double-counting at tile boundaries**
  - **Location:** `src/shared/fragment_iterator.rs:219-223`
  - **Severity:** LOW
  - **Description:** Known issue documented in a TODO: fragment counters can double-count fragments that fall in the halo regions of multiple tiles. The actual analysis output is correct (fragments are correctly assigned to tiles), but the reported statistics (total reads/fragments processed) may be slightly inflated. Either fix or document as a known limitation.

---

- [ ] **T4-11: Improve error messages in reducer**
  - **Location:** `src/commands/fcoverage/reducer.rs` (multiple locations)
  - **Severity:** LOW
  - **Description:** Error messages like `"Missing orig_idx in partials"` don't include context about which tile, chromosome, or line caused the failure. When debugging failed runs on large datasets, this makes root cause analysis difficult. Add tile index and chromosome name to error messages.

---

- [ ] **T4-12: Resolve `gc_bias/config.rs:82` design TODO**
  - **Location:** `src/commands/gc_bias/config.rs:82`
  - **Severity:** LOW
  - **Description:** TODO says: `"Try excluding the first N bases (both ends) from GC fraction calculation to avoid correcting 'biochemical cut bias'."` This is a scientific design decision that could affect correction quality. Decide whether to implement before release or document as a known limitation/future improvement.

---

## Summary Statistics

| Tier | Count | Description |
|------|-------|-------------|
| **T1 (Must Fix)** | 6 | License, metadata, generated code validation, README |
| **T2 (Strongly Recommended)** | 13 | Error handling, dead code, test coverage, clippy |
| **T3 (Should Fix)** | 12 | Unwrap/expect panics, bounds checks, edge cases |
| **T4 (Nice to Have)** | 12 | Code quality, maintainability, documentation |
| **Total** | **43** | |

---

## Test Coverage Per Released Command

| Command | Integration Tests | Unit/Component Tests | Coverage Assessment |
|---------|------------------|---------------------|-------------------|
| `fcoverage` | 7 tests | 31+ tests (coverage, tiling, windows) | **Good** |
| `lengths` | Yes (multi-file) | 25+ tests (counting, tiling) | **Good** |
| `gc-bias` | No integration test | 584 lines across 6 test files | **Moderate** |
| `ref-gc-bias` | No integration test | 2 tests (75 lines) | **Weak** |
| `coverage-weights` | **None** | Stride/smoothing unit tests only | **Weak** |
| `midpoints` | **None** | Partial (1 heatmap test, some profile group tests) | **Weak** |
| `bam-to-bam` | 7 tests | None | **Moderate** |
| `bam-to-frag` | 1 smoke test | None | **Weak** |
| `frag-to-bam` | Multi-test file | Unit tests | **Moderate** |

---

## Appendix: Full TODO/FIXME Inventory in Released Command Paths

| File | Line | Comment |
|------|------|---------|
| `shared/frag_file.rs` | 1 | `TODO: Fix this. Just generated...` |
| `shared/fragment_iterator.rs` | 219 | `TODO: ...might end up counting fragments...in multiple tiles` |
| `shared/blacklist/load.rs` | 18 | `TODO: plumb through a streaming reader` |
| `shared/read.rs` | 75 | `TODO: Requires testing!` (stale — tests exist) |
| `shared/bed.rs` | 1081 | `TODO: Generalize and test` |
| `shared/base.rs` | 32 | `TODO: What is "anything else" possibly?` |
| `shared/tiled_run.rs` | 111 | `TODO: Already initialized with None?` |
| `shared/fragment/with_records_fragment.rs` | 98 | `TODO: Avoid cloning` |
| `shared/kmers/kmer_codec.rs` | 130 | `TODO: Calculate actual limit possible!` |
| `commands/gc_bias/interpolation.rs` | 1 | `TODO: Validate that it's correct` |
| `commands/gc_bias/config.rs` | 82 | `TODO: Try excluding the first N bases...` |
| `commands/gc_bias/gc_bias.rs` | 103 | `TODO: Rename to something meaningful` |
| `commands/gc_bias/gc_bias.rs` | 355 | `keep_temp = false` (dead code) |
| `commands/gc_bias/gc_bias.rs` | 451 | `TODO: not currently used downstream` |
| `commands/gc_bias/gc_bias.rs` | 579 | `TODO: Update this pipeline list` |
| `commands/ref_gc_bias/config.rs` | 4 | `TODO: Do we need to add end-offset here...` |
| `commands/coverage_weights/config.rs` | 5 | `TODO: Improve docstring` |
| `commands/coverage_weights/striding.rs` | 24 | `TODO: Should not happen?` |
| `commands/fcoverage/tiling.rs` | 193 | `TODO: Rename variables to meaningful names` |
| `commands/fcoverage/config.rs` | 195 | `TODO: Consider whether blacklist is "filtering"...` |
| `commands/fcoverage/fcoverage.rs` | 523 | `keep_temp = false` (dead code) |
| `commands/lengths/lengths.rs` | 355 | `keep_temp = false` (dead code) |
| `commands/midpoints/counting_by_group.rs` | 335 | `TODO: Test!!` |
| `commands/midpoints/midpoints.rs` | 267 | `keep_temp = false` (dead code) |
| `commands/bam_to_frag/bam_to_frag.rs` | 259 | `TODO: Consider tiling...to decrease memory` |
| `commands/cli_common.rs` | 274 | `TODO: ...add window-based overlap variants` |
| `commands/cli_common.rs` | 358 | `TODO: Standardize whether lists should be comma-sep or space-sep` |
| `commands/cli_common.rs` | 549 | `TODO: Is "nearest" clear enough...` |
