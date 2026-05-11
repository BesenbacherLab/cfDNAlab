# cfDNAlab Release Readiness Report

**Date:** 2026-04-19
**Version:** 0.1.0
**Branch:** `main`
**Target:** Public GitHub release, crates.io publication, scientific conference + LinkedIn announcement

---

## Executive Summary

| Area               | Status  | Notes                                                                    |
| ------------------ | ------- | ------------------------------------------------------------------------ |
| **Build**          | PASS    | Clean compile (debug + release), 0 errors                                |
| **Tests**          | PASS    | 1 080+ tests, all passing                                                |
| **Clippy**         | PASS    | 0 clippy warnings                                                        |
| **License**        | PASS    | MIT license file present, `license = "MIT"` in Cargo.toml                |
| **Cargo metadata** | PARTIAL | Keywords exceed crates.io limit (8/5)                                    |
| **README**         | BLOCKED | 6 visible TODO placeholders remain                                       |
| **Generated code** | BLOCKED | 2 modules carry "unvalidated AI-generated" warnings                      |
| **Website docs**   | PARTIAL | Missing guides for `ends`, `fragment-count-weights`, conversion commands |

**Overall verdict:** Close to release but with a handful of blockers and several should-fix items.

---

## BLOCKERS

These must be resolved before public release. Each one is visible to users, reviewers, or the crates.io registry and would undermine trust.

### B1. README contains 6 visible TODO placeholders

Users landing on the GitHub page will see unfinished sentences:

| Line | Content                                             |
| ---- | --------------------------------------------------- |
| 9    | "The package is in alpha-stage (being developed)"   |
| 75   | `[TODO: Not that simple]` and `(TODO on samtools!)` |
| 343  | `[TODO: Note on how to get griffin-like profiles]`  |
| 351  | `[TODO: Add example]` (end motifs recipe)           |
| 355  | `[TODO: Add output-prefix for remaining commands]`  |
| 445  | `[TODO: Check correctness]` for file column docs    |

**Action:** Fill in or remove every TODO. Change line 9 to release-ready language. Add the missing end-motifs recipe or remove the placeholder section.

### B2. Source module carries "unvalidated AI-generated code" warning

| File                                      | Warning                              |
| ----------------------------------------- | ------------------------------------ |
| `src/commands/gc_bias/interpolation.rs:1` | `"TODO: Validate that it's correct"` |
| ~~`src/shared/frag_file.rs:1`~~           | ~~Removed~~                          |

The `interpolation.rs` comment will be visible to anyone inspecting the source and is damaging for a tool targeting scientific credibility. Even if the code is now correct, the warning must be removed or replaced with a validation note before release.

**Action:** Validate `interpolation.rs` (or confirm it's already been validated) and remove or rewrite the warning comment.

### B3. `keywords` in Cargo.toml exceeds the crates.io limit

crates.io allows a maximum of **5 keywords**. The current list has **8**:
`["cfdna", "cell-free-dna", "fragmentomics", "fragmentation", "bioinformatics", "genomics", "sequencing", "cli"]`

`cargo publish` will fail.

**Action:** Reduce to 5 keywords. Suggested: `["cfdna", "cell-free-dna", "fragmentomics", "bioinformatics", "genomics"]`.

---

## SHOULD FIX

Not strict blockers, but each one poses a risk to credibility, usability, or CI reliability for a public release.

### S2. CHANGELOG needs updating

The CHANGELOG currently:
- Sources a `release-notes.md` that says "Multiple other commands are being built. While they are technically present... they either not tested properly or outright won't work." This reads poorly for a release announcement.
- Lists `cfdna ends` as a public command — confirm this is intentional for 0.1.0.
- Lists `cfdna fragment-count-weights` — this is a new command not present in the older reviews. Confirm it's release-ready.

**Action:** Rewrite the CHANGELOG release notes for 0.1.0 to be professional and confident. Remove the "won't work" language.

### S3. Website intro doesn't mention `ends` or `fragment-count-weights`

The [intro.md](website/docs/intro.md) says cfDNAlab extracts "fragment coverage, midpoint coverage, and fragment lengths" — it omits **fragment end- and breakpoint motifs**.

**Action:** Add end-motifs to the intro feature list, matching the README.

### S5. `delfi_features_guide.md` has 5 unresolved TODOs

This is the only website guide with TODOs. All five are research questions about how DELFI actually works. If this guide ships as-is, users will see incomplete documentation.

**Action:** Either finish the guide or remove/hide it from the public docs until it's complete.

### S6. Dead `keep_temp = false` pattern in 5 commands

Five commands contain `let keep_temp = false;` with a dead `else` branch:
- `lengths.rs:467`, `gc_bias.rs:397`, `midpoints.rs:302`, `ends.rs:431`, `wps.rs:607`

This is not user-facing but looks sloppy in source review. 

**Action:** Remove the dead branches or implement `--keep-temp` behind a dev feature.

### ~~S6. Commented-out code~~ — RESOLVED

Commented-out functions in `bam.rs` and `striding.rs` removed by user. `frag_file.rs` warning removed. `lib.rs` re-exports are intentional (future API curation).

---

## NICE TO HAVE

These improve polish but won't block a successful release.

### N1. Source-level TODOs in released command paths

48 `TODO`/`FIXME` comments remain across 28 source files. Most are minor (naming suggestions, possible optimizations, speculative features). None are user-facing. The count is down from ~65 in March, so progress is real.

Not a blocker — internal TODOs are normal in actively developed software. However, the two AI-generated-code TODOs (B2) stand out and must be addressed.

### N2. Missing website guides

Guides exist for: fragment coverage, fragment lengths, fragment midpoints, genomic smoothing, GC bias correction, DELFI features.

Guides missing for: **end motifs**, **fragment-count-weights**, **bam-to-bam**, **bam-to-frag**, **frag-to-bam** (conversion commands).

The README recipes partially cover some of these. For a first release, the README recipes + CLI `--help` may be sufficient, but guides would strengthen the documentation story for a conference presentation.

### N3. Missing CLI reference docs for `ends` and `fragment-count-weights`

The generated CLI docs in `website/docs/generated/cli/` cover 10 commands but are missing `ends` and `fragment-count-weights`. The doc generation script likely needs updating.

---

## WHAT'S IN GOOD SHAPE

These areas have been resolved since the earlier reviews and are now release-ready:

1. **Clippy:** 0 warnings (was 425 in March). This is excellent for a public release.

2. **Cargo.toml metadata:** `authors`, `description`, `license`, `repository`, `homepage`, `documentation`, `categories` are all populated and correct.

3. **LICENSE file:** MIT license present with correct copyright.

4. **CI infrastructure:** Three GitHub Actions workflows exist:
   - `rust.yml` — build + test with `--all-features`
   - `code_cov.yml` — Codecov coverage with `--all-features`
   - `docs.yml` — generates CLI docs, builds and deploys Docusaurus site to GitHub Pages

5. **Test coverage:** 1 080+ test functions across integration and unit tests. All commands have substantial test coverage. The test plan identifies no high-priority gaps remaining.

6. **Test quality:** The test suite covers fragment semantics, window semantics, tiling invariance, weight composition, cross-command interoperability, and artifact contracts. A release-spine workflow test verifies the full `ref-gc-bias -> gc-bias -> coverage-weights -> fcoverage` pipeline.

7. **CONTRIBUTING.md:** Present at `.github/CONTRIBUTING.md`.

8. **Website documentation:** Docusaurus site with auto-generated CLI reference, user guides, and installation instructions. Deployed via GitHub Pages.

9. **Tile/window refactor:** The April 2 refactor (BedFetchPolicy, candidate window spans) is architecturally sound. Boundary logic has been manually verified. The main risk (fcoverage double-counting from wrong model) is now tested.

---

## RELEASE CHECKLIST

### Before tagging 0.1.0

- [ ] Remove all README TODOs (B1)
- [ ] Remove/rewrite AI-generated-code warnings in `interpolation.rs` (B2)
- [X] Reduce Cargo.toml keywords to 5 (B3)
- [X] Rewrite CHANGELOG to be release-quality (S2)
- [X] Update website intro (S3)
- [ ] Fix or remove DELFI guide TODOs (S4)
- [ ] Remove dead `keep_temp` branches (S5)
- [X] Remove commented-out code (S6)

### Before `cargo publish`

- [ ] Run `cargo package --features cli` and verify it succeeds
- [ ] Verify `cargo publish --dry-run` passes
- [ ] Confirm the binary name `cfdna` doesn't conflict on crates.io

### Before conference/LinkedIn announcement

- [ ] Verify website deploys correctly (https://BesenbacherLab.github.io/cfdnalab/)
- [ ] Add end-motifs guide or at minimum a complete README recipe (N2)
- [ ] Update generated CLI docs to include `ends` and `fragment-count-weights` (N3)
- [ ] Prepare a brief feature summary suitable for LinkedIn (command list, speed claims, key differentiators)

---

## CLI DOCS REVIEW

Reviewed `--help` output for all 12 commands listed in the README plus the top-level `cfdna --help`. Findings are organized by severity.

### Inconsistencies across commands

| #   | Finding                                                                                                                                                                                       | Commands                                     | Severity |
| --- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------- | -------- |
| C8  | **`--output-prefix` examples inconsistent.** `coverage-weights` and `fragment-count-weights` show no example output filenames. Other commands do.                                             | `coverage-weights`, `fragment-count-weights` | LOW      |
| C9  | **`--require-proper-pair` wording differs.** `fcoverage` uses a shorter version without the "trims the tails" explanation. Others include it. Minor, but noticeable when comparing help text. | `fcoverage` vs others                        | LOW      |

---

## DOCUMENTS SUPERSEDED BY THIS REPORT

The following documents are consolidated into this report and should be removed:

| File                                                       | Original date | Reason for removal                                                                      |
| ---------------------------------------------------------- | ------------- | --------------------------------------------------------------------------------------- |
| `CURRENT_STATE_REVIEW.md`                                  | 2026-04-02    | Tile/fetch refactor review — findings incorporated, code has evolved                    |
| `RELEASE_READINESS_REPORT.md`                              | 2026-03-12    | Original readiness report — most items resolved or superseded                           |
| `RELEASED_COMMANDS_TEST_PLAN.md`                           | 2026-03-25    | Test plan — coverage gaps closed, plan self-reports "no concrete high-priority gaps"    |
| `plans_and_specs/RELEASE_TODO.md`                          | Various       | Release TODO — blockers tracked here instead                                            |
| `plans_and_specs/COMMAND_RELEASE_AUDIT_SPEC_2026-03-03.md` | 2026-03-03    | Audit spec — findings resolved or tracked here                                          |
| `plans_and_specs/ends_code_review.md`                      | Various       | Ends code review — `ends` is now substantially tested; remaining items are low-severity |
