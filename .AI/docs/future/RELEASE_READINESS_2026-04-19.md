# cfDNAlab Release Readiness Review

**Date:** 2026-05-15
**Version:** 0.1.0
**Target:** first public release, crates.io package, public docs, announcement

This replaces the older 2026-04-19 review. The old README TODO, Cargo keyword, generated-doc, and AI-warning blockers are no longer the right release focus.

## Verdict

The release now looks like a final-verification problem, not a broad cleanup problem. The package boundary, public docs URL, and DELFI draft exposure have been addressed in the repo. Do not spend release-day time on wider refactors.

## Necessary Before Release

### 1. Verify the crates.io package artifact

`Cargo.toml` has an explicit package whitelist:

- `/Cargo.toml`
- `/LICENSE`
- `/README.md`
- `/cfdnalab_logo_750x500_150dpi.png`
- `/src/**/*.rs`
- `/tests/**/*.rs`

This keeps the crate focused on code, the main README, the README image, and the test suite. It excludes `.AI/`, `website/`, `.github/`, `other_software/`, generated site output, and Node dependencies.

Action:

- Run `cargo package --list` after committing, or `cargo package --list --allow-dirty` before committing, and check that no internal/docs/vendor paths appear.
- Run `cargo package`.
- Run `cargo publish --dry-run`.

### 2. Verify custom-domain docs deployment

The canonical public docs URL is `https://cfDNAlab.tools/`.

Current repo settings point there:

- `README.md` uses `https://cfDNAlab.tools` / `https://cfDNAlab.tools/` for docs links.
- `Cargo.toml` uses `https://cfDNAlab.tools/` for homepage and documentation.
- `website/docusaurus.config.js` uses `url: 'https://cfDNAlab.tools'` and `baseUrl: '/'`.

Action:

- After deployment, verify `https://cfDNAlab.tools/` loads the built site.
- Verify docs links resolve from the site root.
- Verify the GitHub Pages custom domain and DNS settings are configured in GitHub/DNS. With the current GitHub Actions Pages workflow, a `CNAME` file is not required.

### 3. Run the final checks

I did not run tests. Earlier compile-only checks completed cleanly during this review, but you should rerun final checks yourself from the release state:

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

Then run your normal full test suite.

## Maybe Necessary

### Cargo.lock is ignored and untracked

`Cargo.lock` exists locally but is ignored by `.gitignore` and not tracked. For a CLI distributed via `cargo install --git` and crates.io, this is a reproducibility decision.

Action:

- Decide whether the release should commit `Cargo.lock`.
- If yes, remove `Cargo.lock` from `.gitignore`, commit the lockfile, and use locked installs/checks where relevant.
- If no, accept that installs and CI resolve current compatible dependency versions.

### README badges need post-publish cleanup

The README badges still say crates.io and BioConda are not published. This is acceptable for the initial release, but keep it in the follow-up plan.

Action:

- After crates.io publication, update the crates.io badge.
- Do not claim BioConda availability until it exists.
- Revisit badges for the GitHub repo and v0.2 cleanup.

### CHANGELOG is minimal

`CHANGELOG` is release-safe now. It lists the public commands and no longer contains the old "won't work" language. It is still thin for an announcement.

Action:

- Optional: add a short 0.1.0 summary with the command groups, supported input model, and bias/smoothing support.
- Not necessary if the README or announcement text carries that context.

## Accepted Risks

- `.github/workflows/rust_all_features.yml` intentionally signals failure if the all-features build/test job fails. The comment about not interpreting every failure as a release blocker is not a release issue.

## No Longer Relevant From The Old Report

- README TODO placeholders: not present in the current README.
- Cargo keyword limit: fixed at five keywords.
- `src/commands/gc_bias/interpolation.rs` AI-validation warning: removed.
- Missing generated CLI pages for `ends` and `fragment-count-weights`: current generated docs include both.
- Dead `keep_temp = false` branches: no current matches in `src/commands` or `src/shared`.
- Website intro missing `ends`: fixed.
- DELFI guide exposure: moved out of `website/docs/` to `website/drafts/delfi_features_guide.md`.
- Source-visible generated/unvalidated test TODOs: removed from `tests/test_sampling.rs` and `tests/test_prepare_windows_near.rs`.

## Current Good Signals

- `Cargo.toml` metadata has name, version, authors, description, license, README, repository, homepage, documentation, keywords, and category.
- Default features expose the documented release commands and do not include experimental commands such as `fragment-kmers`, `prep-windows`, `wps`, or `wps-peaks`.
- README command list matches the default release command set.
- Generated CLI docs exist locally for the 11 release commands and the generation script uses `--scope release`.
- CI has separate default, all-features, docs, and coverage workflows.
- The Docusaurus config uses the custom domain and root base URL.

## Release-Day Order

1. Verify the package whitelist with `cargo package --list`.
2. Run `cargo package` and `cargo publish --dry-run`.
3. Run the final compile/docs checks listed above.
4. Run your full tests.
5. Verify `https://cfDNAlab.tools/` after docs deployment.
6. Publish/tag.
7. Update badges and announcement links after packages/docs are live.
