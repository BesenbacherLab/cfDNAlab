# website Spec

The website is a Docusaurus documentation site whose command reference is generated from the Rust Clap command tree.

## Site Stack

- Website root is `website/`.
- Docusaurus version is pinned through `website/package.json`.
- Node must be at least 20.
- Build output goes to `website/generated-site`.
- `npm run start` and `npm run build` regenerate CLI docs first through package pre-scripts.

## Generated CLI Reference

- Generator binary is `gen_cli_docs` and requires features `cli,docs_gen`.
- Generated CLI pages live under `website/docs/generated/cli/`.
- Generated files include marker comments:

```text
<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->
```

- `GENERATED_NOTICE.txt` records the regeneration command.
- `--scope release` includes only release commands:
  - `bam-to-bam`
  - `bam-to-frag`
  - `coverage-weights`
  - `ends`
  - `fcoverage`
  - `fragment-count-weights`
  - `frag-to-bam`
  - `gc-bias`
  - `lengths`
  - `midpoints`
  - `ref-gc-bias`
- `--scope all` includes every command exposed by the built command tree.

## Generator Behavior

- The generator removes stale generated Markdown pages before writing new ones.
- Help text comes from Clap long help.
- Host-dependent `--n-threads` defaults are normalized to `[default: auto]`.
- Generated overview and command pages are sorted by command name.
- The parser keeps long-about prose as intro markdown, extracts usage, and groups options by help heading.
- Nested bullet prose is indented so it renders under the relevant option.
- `--fail-on-drift` runs `git diff --exit-code` on the generated output directory.

## Release Notes

- `website/scripts/generate_release_notes.sh` writes `website/docs/generated/release-notes.md`.
- Source is `CHANGELOG`.
- The generated page has an auto-generated marker and removes the top-level changelog title.

## CI Contract

- `.github/workflows/docs.yml` builds docs on pull requests and pushes to `main`.
- The docs build job regenerates CLI docs and release notes before installing Node dependencies and running `npm run build`.
- Generated docs source files are build inputs in the GitHub Actions runner and are not required to be committed.
- Pull requests validate the docs build without uploading a site artifact.
- Pushes to `main` upload the built site artifact and deploy it to GitHub Pages.

## Editing Rules

- Do not hand-edit files under `website/docs/generated/cli/` or generated release notes.
- Update Rust Clap config and long help when command behavior changes, then regenerate docs.
- Update authored guide pages under `website/docs/guides/` when behavior needs narrative explanation beyond CLI reference.
- Keep generated-doc feature lists in scripts and `GENERATED_NOTICE.txt` in sync with release command features.
