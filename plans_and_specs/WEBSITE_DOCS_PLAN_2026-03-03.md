# cfDNAlab Docs Website Plan

Date: 2026-03-03

## Summary

Build a Docusaurus website in this repo and deploy it to GitHub Pages from GitHub Actions.

CLI help pages must be auto-generated from the same Clap command definitions used by `cfdna`.

All generated content must be placed in clearly marked generated folders, with CI checks that prevent manual edits from drifting.

## Locked decisions

1. Hosting platform: GitHub Pages
2. Docs scope at launch: 9 release commands only
3. Versioning at launch: latest only
4. Website framework: Docusaurus
5. Generated-content policy:
- Generated CLI markdown pages are committed for review and diff visibility
- Static site build output is not committed
6. Tooling isolation policy:
- Website tooling is isolated in `website/` and must not be required for Rust build/test flows
- Any new Rust dependencies for docs generation must be optional and gated behind a non-default feature

## Goals and non-goals

### Goals

1. Pretty docs site with search
2. Custom front page with logo, package intro, and clear links
3. Auto-generated CLI reference pages from Clap definitions
4. CI-driven regeneration and deploy on `main`
5. Strong generated-folder safeguards

### Non-goals for v1

1. Multi-version docs
2. Algolia integration
3. Publishing experimental commands by default

## Source-of-truth model

1. CLI docs source of truth:
- Rust Clap docstrings in command config structs and command tree definitions
2. Website pages source of truth:
- Manual docs under `website/docs/` (non-generated)
- Generated CLI pages under `website/docs/generated/cli/`

No manual editing is allowed in generated folders.

## Repository layout

Create and maintain the following structure:

1. `website/`
- Docusaurus project root

2. `website/docs/`
- Manual docs

3. `website/docs/generated/cli/`
- Auto-generated command pages only
- Must include warning marker in each file

4. `website/static/img/`
- Logo and static assets

5. `website/.generated-site/`
- Docusaurus build output folder
- Gitignored
- Generated only by build

6. `scripts/docs/`
- Small helper scripts for docs generation checks

7. `src/cli/` (new shared module)
- Shared Clap command-building functions used by both binary and docs generator

8. `src/bin/gen_cli_docs.rs`
- CLI docs generator binary

9. No root-level Node tooling
- No repo-root `package.json`, lockfile, or Node scripts
- Node dependencies live only in `website/`

## Interface and API changes

## Dependency isolation rules

1. Rust dependencies used only by docs generation must be declared `optional = true`.
2. Add a dedicated non-default feature, for example `docs_gen`, to gate those dependencies.
3. `gen_cli_docs` binary must require docs features explicitly:
- `required-features = ["cli", "docs_gen"]`
4. `default` Cargo features must remain unchanged and must not include `docs_gen`.
5. CI docs job is the only required path that enables `docs_gen`.

## 1) Shared CLI builder in library

Add shared command-building APIs so docs generator and runtime use the same command graph:

1. `pub fn build_cli_command_for_terminal() -> clap::Command`
2. `pub fn build_cli_command_for_docs() -> clap::Command`

Rationale:
- Terminal builder keeps current presentation behavior
- Docs builder avoids terminal decoration artifacts that do not belong in markdown docs

## 2) Generator binary interface

Create `gen_cli_docs` binary with this contract:

Command:

```bash
cargo run --bin gen_cli_docs --features cli,cmd_bam_to_bam,cmd_bam_to_frag,cmd_frag_to_bam,cmd_coverage_weights,cmd_fcoverage,cmd_gc_bias,cmd_lengths,cmd_midpoints,cmd_ref_gc_bias -- --out-dir website/docs/generated/cli --scope release
```

If docs-only dependencies are added, include `docs_gen` in that invocation.

Supported flags:

1. `--out-dir <path>`
2. `--scope <release|all>` (default `release`)
3. `--fail-on-drift` (optional CI validation mode)

## 3) Generated file contract

Every generated markdown file must begin with:

```md
<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->
```

`website/docs/generated/cli/README.md` must also explain:

1. Why files are generated
2. Exact regeneration command
3. No-manual-edit policy

## Website content and UX spec

## 1) Front page

Custom homepage (`website/src/pages/index.tsx`) includes:

1. Logo from `website/static/img/`
2. One-paragraph package description
3. Primary actions:
- Get started
- Commands
- Release notes
4. Clear top nav and section hierarchy

Theme and styling:

1. Custom CSS in `website/src/css/custom.css`
2. Branded color palette
3. Mobile-first layout checks

## 2) CLI reference section

Add sidebar section `CLI Reference` with one generated page per release command:

1. `gc-bias`
2. `ref-gc-bias`
3. `coverage-weights`
4. `lengths`
5. `fcoverage`
6. `midpoints`
7. `bam-to-bam`
8. `bam-to-frag`
9. `frag-to-bam`

Also generate an index page in `website/docs/generated/cli/index.md`.

## 3) Search

Use Docusaurus local search plugin at launch.

Keep structure compatible with future Algolia adoption without content migration.

## Generation behavior details

Generation pipeline rules:

1. Use `build_cli_command_for_docs()` as the command tree source
2. Traverse command and subcommand structure deterministically
3. Render markdown with deterministic ordering of sections and options
4. Write only to `website/docs/generated/cli/`
5. Never write into non-generated docs paths
6. Produce stable output across repeated runs

## No-manual-edit safeguards

Implement all safeguards below:

1. Folder naming:
- Generated docs must stay under `website/docs/generated/`

2. In-file markers:
- Every generated file includes no-edit marker header

3. CI drift check:
- Regenerate docs in CI
- Fail when `git diff --exit-code website/docs/generated/cli` is non-zero

4. Contributor guidance:
- Add short rule to contributor docs stating generated paths are not edited manually

5. Build output isolation:
- Docusaurus output path set to `website/.generated-site/` and gitignored

## GitHub Actions workflow

Add `.github/workflows/docs.yml` with three jobs:

## 1) `generate-cli-docs`

1. Checkout
2. Install Rust toolchain
3. Run generator with release feature set and `docs_gen` feature
4. Verify markers exist in generated files
5. Fail if generated files are stale (`git diff --exit-code`)

## 2) `build-site`

1. Setup Node
2. Install dependencies in `website/`
3. Build Docusaurus site to `website/.generated-site/`
4. Upload artifact for PR visibility

## 3) `deploy-pages` (push to `main` only)

1. Use GitHub Pages official actions
2. Deploy built artifact

Triggers:

1. Pull request to `main`: generate + build (no deploy)
2. Push to `main`: generate + build + deploy

## Implementation phases

## Phase 1: CLI sharing and generator foundation

1. Extract shared command-tree builder from current binary path
2. Implement `gen_cli_docs` binary
3. Generate release command pages to `website/docs/generated/cli/`
4. Add file markers and generated folder README

## Phase 2: Docusaurus scaffold and front page

1. Initialize `website/` project
2. Add homepage with logo, description, and links
3. Add docs navigation and CLI reference wiring
4. Add custom theme styling

## Phase 3: CI and deployment

1. Add docs GitHub Actions workflow
2. Add drift checks and generated marker checks
3. Deploy to GitHub Pages on `main`

## Phase 4: hardening

1. Improve command page readability template
2. Validate local search behavior
3. Add contributor docs for docs workflow

## Test and validation plan

## Generator tests

1. Generates all 9 release command pages
2. Includes required marker header in each generated file
3. Output is deterministic across repeated runs
4. Includes expected sections such as usage and options

## Integration checks

1. Editing a command docstring changes generated page diff
2. CI fails when generated pages are stale
3. CI passes when generation is up to date

## Site checks

1. Homepage renders with logo and key links
2. CLI pages are reachable from sidebar
3. Search returns command pages
4. Mobile layout is usable

## Acceptance criteria

1. Docs site is live on GitHub Pages
2. Front page is customized and branded
3. CLI reference pages for all 9 release commands are generated automatically
4. Generated folders are clearly marked and protected by CI
5. Changing command help in Rust updates website content via CI regeneration and deploy

## Assumptions and defaults

1. Release commands remain the default docs scope
2. Experimental commands are excluded at launch
3. Latest-only docs at launch
4. Generated markdown is committed, static built site is not
5. Clap docstrings remain the single source of truth for CLI docs
