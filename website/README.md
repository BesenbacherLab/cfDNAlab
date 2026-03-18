# cfDNAlab docs website

This folder contains the Docusaurus website for cfDNAlab documentation.

## Local development

1. Install dependencies:

```bash
npm install
```

2. Start local dev server:

```bash
npm run start
```

`npm run start` generates the CLI reference docs and release notes before starting Docusaurus.

Open the URL printed by Docusaurus in your browser.

## Build static site

```bash
npm run build
```

`npm run build` generates the CLI reference docs and release notes before building the site.

Build output is written to `website/generated-site/`.

## Generated folders

- `website/docs/generated/cli/` is auto-generated and must not be edited manually.
- `website/docs/generated/release-notes.md` is auto-generated from `CHANGELOG`.
- `website/generated-site/` is generated build output and must not be committed.
