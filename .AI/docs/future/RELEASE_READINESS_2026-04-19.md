# cfDNAlab Release Readiness Review

**Date:** 2026-05-15
**Version:** 0.1.0
**Target:** first public release, crates.io package, public docs, announcement

This replaces the older 2026-04-19 review. The old README TODO, Cargo keyword, generated-doc, and AI-warning blockers are no longer the right release focus.

## Verdict

The source tree is close, but I would not publish the crate or announce the docs until the package contents, incomplete docs page, and public URL story are fixed or explicitly accepted.

For a same-day release, do not spend time on broad cleanup. Fix the necessary items below, run the release checks locally, then tag.

## Necessary Before Release

### 1. Verify the crates.io package contents

`Cargo.toml` now has an explicit `include` whitelist for the crate package:

- `/Cargo.toml`
- `/LICENSE`
- `/README.md`
- `/cfdnalab_logo_750x500_150dpi.png`
- `/src/**/*.rs`

This should keep the package focused on the CLI code, library code, license, main README, and the image referenced by the README.

The reason this matters is that the repo has tracked internal material that should not go to crates.io by accident:

- `.AI/` contains 50 tracked files, including future plans, major reviews, and this readiness report.
- `website/` contains authored docs and Docusaurus build inputs.
- `.github/` contains CI files that are useful in the repo but not needed in the crate archive.

Action:

- Verify with `cargo package --list --allow-dirty` before committing, or `cargo package --list` after committing, that `.AI/`, old review notes, future plans, website docs, and `.github/` are not included.
- Run `cargo package` after the package list looks right.

This is the highest-impact release hygiene issue because it affects the artifact that gets permanently distributed.

### 2. Remove or hide the unfinished DELFI guide

`website/docs/guides/delfi_features_guide.md` is still under the Docusaurus docs tree and contains visible TODOs plus unfinished prose. It is not in `sidebars.js`, but it is still a docs file and the site config indexes docs content.

Current visible problems include:

- `[TODO: Make it clear that we don't do exactly what they do ...]`
- `[TODO: Figure out exactly what the DELFI settings are]`
- `[TODO: The actual DELFI approach needs to be found in their code]`
- The file ends mid-thought after "Since we already have the fragment-level GC correction, we will".

Action:

- Move this guide out of `website/docs/` until it is complete, or mark it as non-published if the current Docusaurus setup supports that cleanly.
- Do not spend release-day time finishing the DELFI method unless it is central to the announcement.

### 3. Decide the canonical docs URL

The repo currently points users to two different docs homes:

- `README.md` badge and prose use `https://cfDNAlab.tools/`.
- `Cargo.toml` uses `https://BesenbacherLab.github.io/cfdnalab/` and `/docs`.
- `website/docusaurus.config.js` is configured for GitHub Pages with `url: 'https://BesenbacherLab.github.io'` and `baseUrl: '/cfdnalab/'`.

Action:

- If `cfDNAlab.tools` is the public URL, configure the site and Cargo metadata for that URL and make sure deployment writes the right domain setup.
- If GitHub Pages is the public URL, update the README badge and docs links.

Do this before crates.io publication and announcement, because crate metadata and public links are what new users will copy first.

### 4. Run the release checks locally

I started compile/package checks during this review before you redirected the scope. The compile-only checks that completed were clean:

- `cargo check`
- `cargo check --all-features`
- `cargo check --tests`

I did not run tests. `cargo package` could not complete in the sandbox because it needed crates.io index access, and the later offline attempt was interrupted.

Action for you before tagging/publishing:

```bash
cargo check
cargo check --all-features
cargo check --tests
cargo package --list
cargo package
cargo publish --dry-run
npm --prefix website ci
npm --prefix website run build
```

Then run your normal test suite outside this agent session.

## Maybe Necessary

### Cargo.lock is ignored and untracked

`Cargo.lock` exists locally but is ignored by `.gitignore` and not tracked. For a CLI distributed via `cargo install --git` and crates.io, this is a reproducibility decision, not a style issue.

Action:

- Decide whether the release should commit `Cargo.lock`.
- If yes, remove `Cargo.lock` from `.gitignore`, commit the lockfile, and use locked installs/checks where relevant.
- If no, accept that installs and CI resolve current compatible dependency versions.

### README badges need post-publish cleanup

The README badges still say crates.io and BioConda are not published. That is fine before publication, but it will look stale immediately after the release if the announcement sends users to GitHub.

Action:

- Update package badges after each package is actually published.
- Do not claim BioConda availability until it exists.

### CHANGELOG is minimal

`CHANGELOG` is release-safe now. It lists the public commands and no longer contains the old "won't work" language. It is still thin for an announcement.

Action:

- Optional: add a short 0.1.0 summary with the command groups, supported input model, and bias/smoothing support.
- Not necessary if the README or announcement text carries that context.

### Source-visible unvalidated-test comments remain

There are still source-visible comments in tests such as:

- `tests/test_sampling.rs`: "TODO: Validate tests - generated but not yet checked!"
- `tests/test_prepare_windows_near.rs`: "TODO: These tests are completely unvalidated."

These are not runtime blockers and they are outside the main release docs, but they are visible in a public source review.

Action:

- Leave them if time is tight.
- Clean them if the release goal includes source credibility for scientific reviewers.

### Stale all-features CI comment

`.github/workflows/rust_all_features.yml` says "This is allowed to fail, as not all commands are done yet", but the workflow does not actually allow failure. That comment reads stale for a release tree.

Action:

- Remove or update the comment when convenient.
- This is not a release blocker if the workflow passes.

## No Longer Relevant From The Old Report

- README TODO placeholders: not present in the current README.
- Cargo keyword limit: fixed at five keywords.
- `src/commands/gc_bias/interpolation.rs` AI-validation warning: removed. The file now has ordinary implementation docs and interpolation tests.
- Missing generated CLI pages for `ends` and `fragment-count-weights`: current generated docs include both.
- Dead `keep_temp = false` branches: no current matches in `src/commands` or `src/shared`.
- Website intro missing `ends`: fixed.

## Current Good Signals

- `Cargo.toml` metadata has name, version, authors, description, license, README, repository, homepage, documentation, keywords, and category.
- Default features expose the documented release commands and do not include experimental commands such as `fragment-kmers`, `prepare-windows`, `wps`, or `wps-peaks`.
- README command list matches the default release command set.
- Generated CLI docs exist locally for the 11 release commands and the generation script uses `--scope release`.
- CI has separate default, all-features, docs, and coverage workflows.

## Release-Day Order

1. Verify the `include` whitelist with `cargo package --list`.
2. Remove or hide the DELFI guide from public docs.
3. Make `cfDNAlab.tools` vs GitHub Pages consistent across README, Cargo metadata, and Docusaurus.
4. Run the release checks listed above.
5. Run your full tests.
6. Publish/tag.
7. Update badges and announcement links after packages/docs are live.
