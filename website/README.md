# cfDNAlab docs website

This folder contains the Docusaurus website for cfDNAlab documentation.

## Local development

1. Generate CLI docs:

```bash
../scripts/docs/generate_cli_docs.sh
```

2. Install dependencies:

```bash
npm install
```

3. Start local dev server:

```bash
npm run start
```

Open the URL printed by Docusaurus in your browser.

## Build static site

```bash
npm run build
```

Build output is written to `website/.generated-site/`.

## Generated folders

- `website/docs/generated/cli/` is auto-generated and must not be edited manually.
- `website/.generated-site/` is generated build output and must not be committed.
