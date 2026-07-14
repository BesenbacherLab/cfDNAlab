#!/usr/bin/env python3
"""Generate Docusaurus pages for the R and Python output loader APIs."""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
import html
import json
import re
import shutil
import textwrap
from pathlib import Path


GENERATED_MARKER = "<!-- AUTO-GENERATED FILE - DO NOT EDIT -->"
GENERATED_SOURCE = "<!-- Source: py-cfdnalab docstrings and r-cfdnalab roxygen .Rd files -->"


@dataclass(frozen=True)
class MethodDoc:
    name: str
    signature: str
    docstring: str


@dataclass(frozen=True)
class PythonSymbolDoc:
    name: str
    kind: str
    signature: str
    docstring: str
    methods: tuple[MethodDoc, ...]


@dataclass(frozen=True)
class RArgumentDoc:
    name: str
    description: str


@dataclass(frozen=True)
class RTopicDoc:
    name: str
    aliases: tuple[str, ...]
    title: str
    usage: str
    arguments: tuple[RArgumentDoc, ...]
    value: str
    description: str
    details: str
    examples: str


@dataclass(frozen=True)
class LoaderPage:
    slug: str
    title: str
    sidebar_label: str
    cli_command: str
    output_file: str
    summary: str
    python_symbols: tuple[str, ...]
    r_topics: tuple[str, ...]
    python_example: str
    r_example: str
    notes: tuple[str, ...]


COMMON_PYTHON_SYMBOLS: tuple[str, ...] = ()
COMMON_R_TOPICS: tuple[str, ...] = ()


LOADER_PAGES = (
    LoaderPage(
        slug="midpoint-profiles",
        title="Midpoint Profiles",
        sidebar_label="Midpoint Profiles",
        cli_command="cfdna midpoints",
        output_file="<prefix>.midpoint_profiles.zarr",
        summary=(
            "Load midpoint profile Zarr stores and extract count arrays or data frames "
            "by group, fragment length bin, and midpoint position."
        ),
        python_symbols=("read_midpoints", "MidpointProfiles"),
        r_topics=(
            "read_midpoints",
            "group_metadata",
            "length_bins",
            "positions",
            "group_idx",
            "length_bin_idx",
            "profile_array",
            "midpoint_data_frame",
            "midpoint_array",
            "schema_version",
        ),
        python_example=textwrap.dedent(
            """\
            import cfdnalab as cfl

            profiles = cfl.read_midpoints("sample.midpoint_profiles.zarr")

            groups = profiles.group_metadata()
            length_bins = profiles.length_bins()
            positions = profiles.positions()

            profile = profiles.data_frame(groups="LYL1", with_lengths=167)
            """
        ),
        r_example=textwrap.dedent(
            """\
            library(cfdnalab)

            profiles <- read_midpoints("sample.midpoint_profiles.zarr")

            groups <- group_metadata(profiles)
            length_bins <- length_bins(profiles)
            positions <- positions(profiles)

            profile <- midpoint_data_frame(
              profiles,
              groups = "LYL1",
              with_lengths = 167
            )
            """
        ),
        notes=(
            "Select groups or length bins before expanding to a data frame when the store is large.",
            "Python indices are zero-based. R indices are one-based.",
        ),
    ),
    LoaderPage(
        slug="end-motif-counts",
        title="End-Motif Counts",
        sidebar_label="End-Motif Counts",
        cli_command="cfdna ends",
        output_file="<prefix>.end_motifs.zarr",
        summary=(
            "Load dense or sparse end-motif count Zarr stores and extract motif count "
            "tables, dense arrays, or sparse matrices."
        ),
        python_symbols=(
            "read_end_motifs",
            "EndMotifCounts",
            "GlobalEndMotifCounts",
            "WindowedEndMotifCounts",
            "GroupedEndMotifCounts",
        ),
        r_topics=(
            "read_end_motifs",
            "storage_mode",
            "row_mode",
            "motifs",
            "motif_idx",
            "has_motif",
            "window_metadata",
            "group_metadata",
            "group_idx",
            "end_motif_data_frame",
            "dense_counts_matrix",
            "dense_counts_vector",
            "dense_corrected_counts_matrix",
            "sparse_counts_matrix",
            "sparse_corrected_counts_matrix",
            "schema_version",
        ),
        python_example=textwrap.dedent(
            """\
            import cfdnalab as cfl

            ends = cfl.read_end_motifs("sample.end_motifs.zarr")

            storage = ends.storage_mode()
            row_mode = ends.row_mode()
            motifs = ends.motifs_metadata()

            motif_counts = ends.data_frame(motifs="_AA")
            motif_matrix = ends.sparse_counts_matrix(motifs="_AA")
            """
        ),
        r_example=textwrap.dedent(
            """\
            library(cfdnalab)

            ends <- read_end_motifs("sample.end_motifs.zarr")

            storage <- storage_mode(ends)
            row_mode <- row_mode(ends)
            motif_table <- motifs(ends)

            motif_counts <- end_motif_data_frame(ends, motifs = "_AA")
            motif_matrix <- sparse_counts_matrix(ends, motifs = "_AA")
            """
        ),
        notes=(
            "Sparse stores stay sparse unless a dense helper explicitly receives a densify option.",
            "Windowed and grouped outputs can filter rows by `max_blacklisted_fraction`.",
        ),
    ),
    LoaderPage(
        slug="reference-kmer-frequencies",
        title="Reference K-mer Frequencies",
        sidebar_label="Reference K-mers",
        cli_command="cfdna ref-kmers",
        output_file="<prefix>.ref_kmers.zarr",
        summary=(
            "Load reference k-mer frequency Zarr stores and extract frequency "
            "tables, dense arrays, sparse matrices, and reconstructed counts."
        ),
        python_symbols=(
            "read_ref_kmers",
            "RefKmerFrequencies",
            "GlobalRefKmerFrequencies",
            "WindowedRefKmerFrequencies",
            "GroupedRefKmerFrequencies",
        ),
        r_topics=(
            "read_ref_kmers",
            "storage_mode",
            "row_mode",
            "motif_axis_kind",
            "kmer_size",
            "canonical",
            "all_motifs",
            "assign_by",
            "motifs",
            "motif_idx",
            "window_metadata",
            "group_metadata",
            "group_idx",
            "reference_contig_footprint",
            "row_scaling_factors",
            "ref_kmer_data_frame",
            "dense_frequencies_matrix",
            "dense_frequencies_vector",
            "sparse_frequencies_matrix",
            "schema_version",
        ),
        python_example=textwrap.dedent(
            """\
            import cfdnalab as cfl

            ref_kmers = cfl.read_ref_kmers("sample.ref_kmers.zarr")

            motifs = ref_kmers.motifs_metadata()
            scaling = ref_kmers.row_scaling_factors()

            frequencies = ref_kmers.sparse_frequencies_matrix(motifs="ACGT")
            counts = ref_kmers.data_frame(motifs="ACGT")
            """
        ),
        r_example=textwrap.dedent(
            """\
            library(cfdnalab)

            ref_kmers <- read_ref_kmers("sample.ref_kmers.zarr")

            motif_table <- motifs(ref_kmers)
            scaling <- row_scaling_factors(ref_kmers)

            frequencies <- sparse_frequencies_matrix(ref_kmers, motifs = "ACGT")
            counts <- ref_kmer_data_frame(ref_kmers, motifs = "ACGT")
            """
        ),
        notes=(
            "Reference k-mer stores contain frequencies. Count helpers reconstruct counts with `row_scaling_factor`.",
            "Sparse stores stay sparse unless a dense helper explicitly receives an allow-densify option.",
        ),
    ),
    LoaderPage(
        slug="length-counts",
        title="Length Counts",
        sidebar_label="Length Counts",
        cli_command="cfdna lengths",
        output_file="<prefix>.length_counts.tsv.zst",
        summary=(
            "Load fragment length-count TSV outputs and return counts, fractions, or "
            "densities as arrays, matrices, vectors, or data frames."
        ),
        python_symbols=(
            "read_lengths",
            "LengthCounts",
            "GlobalLengthCounts",
            "WindowedLengthCounts",
            "GroupedLengthCounts",
        ),
        r_topics=(
            "read_lengths",
            "length_bins",
            "length_bin_idx",
            "length_counts_matrix",
            "length_counts_vector",
            "length_data_frame",
            "window_metadata",
            "group_metadata",
            "group_idx",
        ),
        python_example=textwrap.dedent(
            """\
            import cfdnalab as cfl

            lengths = cfl.read_lengths("sample.length_counts.tsv.zst")

            bins = lengths.length_bins()
            counts = lengths.counts_array(with_length_range=(100, 221))
            fractions = lengths.data_frame(value="fraction")
            """
        ),
        r_example=textwrap.dedent(
            """\
            library(cfdnalab)

            lengths <- read_lengths("sample.length_counts.tsv.zst")

            bins <- length_bins(lengths)
            counts <- length_counts_matrix(lengths, with_length_range = c(100L, 221L))
            fractions <- length_data_frame(lengths, value = "fraction")
            """
        ),
        notes=(
            "Use `with_lengths`, `with_length_range`, or direct length-bin indices to select bins.",
            "For fractions and densities, `denominator` controls whether totals use all bins or selected bins.",
        ),
    ),
)


def main() -> None:
    arguments = parse_arguments()
    repository_root = Path(__file__).resolve().parents[2]
    python_package_dir = repository_root / "py-cfdnalab" / "src" / "cfdnalab"
    r_package_dir = repository_root / "r-cfdnalab"
    python_public_symbols = set(parse_python_all(python_package_dir / "__init__.py"))
    r_public_topics = parse_r_exports(r_package_dir / "NAMESPACE")
    python_docs = load_python_docs(python_package_dir, python_public_symbols)
    r_docs = load_r_docs(r_package_dir, r_public_topics)

    output_dir = arguments.out_dir
    cleanup_generated_dir(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    validate_loader_config(python_public_symbols, python_docs, r_public_topics, r_docs)
    write_category_file(output_dir)
    write_notice(output_dir)
    write_overview_page(output_dir)
    write_loader_pages(output_dir, python_docs, r_docs)
    write_python_api_page(output_dir, python_docs)
    write_r_api_page(output_dir, r_docs)

    page_count = len(LOADER_PAGES) + 3
    print(f"Generated {page_count} output loader page(s) in {output_dir}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate R/Python loader docs")
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory for generated loader docs",
    )
    return parser.parse_args()


def cleanup_generated_dir(output_dir: Path) -> None:
    if not output_dir.exists():
        return
    for child in output_dir.iterdir():
        if child.is_dir():
            shutil.rmtree(child)
        elif child.suffix in {".md", ".mdx", ".json", ".txt"}:
            child.unlink()


def write_category_file(output_dir: Path) -> None:
    category = {
        "label": "Output Loaders",
        "position": 4,
        "link": {"type": "doc", "id": "generated/loaders/overview"},
    }
    (output_dir / "_category_.json").write_text(
        json.dumps(category, indent=2) + "\n",
        encoding="utf-8",
    )


def write_notice(output_dir: Path) -> None:
    notice = (
        "AUTO-GENERATED DIRECTORY - DO NOT EDIT\n"
        "Source: py-cfdnalab docstrings and r-cfdnalab roxygen .Rd files\n\n"
        "Regenerate with:\n\n"
        "python3 website/scripts/generate_loader_docs.py "
        "--out-dir website/docs/generated/loaders\n"
    )
    (output_dir / "GENERATED_NOTICE.txt").write_text(notice, encoding="utf-8")


def write_overview_page(output_dir: Path) -> None:
    cards = []
    for page in LOADER_PAGES:
        cards.append(
            "\n".join(
                [
                    f'<a className="loader-card" href="./{page.slug}">',
                    f"  <span>{html.escape(page.title)}</span>",
                    f"  <small>{html.escape(page.output_file)}</small>",
                    f"  <p>{html.escape(page.summary)}</p>",
                    "</a>",
                ]
            )
        )

    content = f"""---
title: Output Loaders
sidebar_label: Overview
sidebar_position: 1
---

{GENERATED_MARKER}
{GENERATED_SOURCE}

# Output Loaders

The R and Python packages load files created by the `cfdna` command-line tool. They do not install, run, or wrap the CLI.

Use the CLI reference for command flags. Use these pages after you have output files and want to inspect, reshape, or plot them in R or Python.

<div className="loader-card-grid">
{chr(10).join(cards)}
</div>

## Full API Indexes

- [Python API](./python-api): generated from public Python symbols and docstrings.
- [R API](./r-api): generated from roxygen `.Rd` topics for exported R functions and S3 generics.

## Generation Sources

- Python documentation comes from `py-cfdnalab/src/cfdnalab`.
- R documentation comes from `r-cfdnalab/man`, which is generated by roxygen from `r-cfdnalab/R`.
"""
    (output_dir / "overview.mdx").write_text(content, encoding="utf-8")


def write_loader_pages(
    output_dir: Path,
    python_docs: dict[str, PythonSymbolDoc],
    r_docs: dict[str, RTopicDoc],
) -> None:
    for sidebar_position, page in enumerate(LOADER_PAGES, start=2):
        python_table = loader_python_symbol_table(page, python_docs)
        r_table = loader_r_symbol_table(page, r_docs)
        notes = "\n".join(f"- {note}" for note in page.notes)
        content = f"""---
title: {page.title}
sidebar_label: {page.sidebar_label}
sidebar_position: {sidebar_position}
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

{GENERATED_MARKER}
{GENERATED_SOURCE}

# {page.title}

`{page.cli_command}` writes `{page.output_file}`. {page.summary}

This is loader API documentation, not command-line help. See [`{page.cli_command}`](../cli/{command_doc_slug(page.cli_command)}) for command flags and output generation.

## Quick Use

<Tabs groupId="loader-language">
<TabItem value="python" label="Python">

```python
{page.python_example.strip()}
```

</TabItem>
<TabItem value="r" label="R">

```r
{page.r_example.strip()}
```

</TabItem>
</Tabs>

## Notes

{notes}

## Function Index

<Tabs groupId="loader-language">
<TabItem value="python" label="Python">

{python_table}

Full details: [Python API](./python-api#{python_anchor(page.python_symbols[0])}).

</TabItem>
<TabItem value="r" label="R">

{r_table}

Full details: [R API](./r-api#{r_anchor_for_page(page, page.r_topics[0])}).

</TabItem>
</Tabs>
"""
        (output_dir / f"{page.slug}.mdx").write_text(content, encoding="utf-8")


def write_python_api_page(output_dir: Path, python_docs: dict[str, PythonSymbolDoc]) -> None:
    page_links = api_page_links(COMMON_PYTHON_SYMBOLS)
    sections = []
    if COMMON_PYTHON_SYMBOLS:
        sections.append(render_common_python_section(python_docs))
    for page in LOADER_PAGES:
        sections.append(render_python_section(page, python_docs))

    content = f"""---
title: Python API
sidebar_label: Python API
sidebar_position: {len(LOADER_PAGES) + 2}
---

{GENERATED_MARKER}
{GENERATED_SOURCE}

# Python API

Generated from public symbols and docstrings in `py-cfdnalab/src/cfdnalab`.

## Jump To

{page_links}

{chr(10).join(sections)}
"""
    (output_dir / "python-api.mdx").write_text(content, encoding="utf-8")


def write_r_api_page(output_dir: Path, r_docs: dict[str, RTopicDoc]) -> None:
    page_links = api_page_links(COMMON_R_TOPICS)
    sections = []
    if COMMON_R_TOPICS:
        sections.append(render_common_r_section(r_docs))
    for page in LOADER_PAGES:
        sections.append(render_r_section(page, r_docs))

    content = f"""---
title: R API
sidebar_label: R API
sidebar_position: {len(LOADER_PAGES) + 3}
---

{GENERATED_MARKER}
{GENERATED_SOURCE}

# R API

Generated from roxygen `.Rd` topics in `r-cfdnalab/man`. S3 method signatures are shown so mode-specific arguments are visible, but normal user code should call the exported generic.

## Jump To

{page_links}

{chr(10).join(sections)}
"""
    (output_dir / "r-api.mdx").write_text(content, encoding="utf-8")


def api_page_links(common_symbols: tuple[str, ...]) -> str:
    links = []
    if common_symbols:
        links.append("- [Common](#common)")
    links.extend(
        f"- [{page.title}](#{section_anchor(page.title)})" for page in LOADER_PAGES
    )
    return "\n".join(links)


def loader_python_symbol_table(
    page: LoaderPage,
    python_docs: dict[str, PythonSymbolDoc],
) -> str:
    rows = ["| Symbol | Type | Summary |", "| --- | --- | --- |"]
    for symbol_name in page.python_symbols:
        symbol = python_docs[symbol_name]
        rows.append(
            "| "
            f"[`{symbol.name}`](./python-api#{python_anchor(symbol.name)})"
            f" | {symbol.kind} | {table_cell(first_sentence(symbol.docstring))} |"
        )
    return "\n".join(rows)


def loader_r_symbol_table(page: LoaderPage, r_docs: dict[str, RTopicDoc]) -> str:
    rows = ["| Topic | Summary |", "| --- | --- |"]
    for topic_name in page.r_topics:
        topic = r_docs[topic_name]
        rows.append(
            "| "
            f"[`{topic_name}`](./r-api#{r_anchor_for_page(page, topic_name)})"
            f" | {table_cell(first_sentence(topic.description or topic.title))} |"
        )
    return "\n".join(rows)


def render_python_section(
    page: LoaderPage,
    python_docs: dict[str, PythonSymbolDoc],
) -> str:
    rows = ["| Symbol | Type | Summary |", "| --- | --- | --- |"]
    details = []
    for symbol_name in page.python_symbols:
        symbol = python_docs[symbol_name]
        rows.append(
            "| "
            f"[`{symbol.name}`](#{python_anchor(symbol.name)})"
            f" | {symbol.kind} | {table_cell(first_sentence(symbol.docstring))} |"
        )
        details.append(render_python_symbol(symbol))
    return f"""## {page.title}

{page.summary}

{chr(10).join(rows)}

{chr(10).join(details)}
"""


def render_common_python_section(python_docs: dict[str, PythonSymbolDoc]) -> str:
    rows = ["| Symbol | Type | Summary |", "| --- | --- | --- |"]
    details = []
    for symbol_name in COMMON_PYTHON_SYMBOLS:
        symbol = python_docs[symbol_name]
        rows.append(
            "| "
            f"[`{symbol.name}`](#{python_anchor(symbol.name)})"
            f" | {symbol.kind} | {table_cell(first_sentence(symbol.docstring))} |"
        )
        details.append(render_python_symbol(symbol))
    return f"""## Common

Shared Python loader API entries that do not belong to one output type.

{chr(10).join(rows)}

{chr(10).join(details)}
"""


def render_python_symbol(symbol: PythonSymbolDoc) -> str:
    anchor = python_anchor(symbol.name)
    if symbol.kind == "function":
        return f"""### `{symbol.name}` {{#{anchor}}}

{render_python_signature_block(symbol.signature)}

{docstring_to_markdown(symbol.docstring)}
"""

    method_index = ""
    method_details = ""
    if symbol.methods:
        rows = ["| Method | Summary |", "| --- | --- |"]
        detail_blocks = []
        for method in symbol.methods:
            method_anchor = python_anchor(f"{symbol.name}.{method.name}")
            rows.append(
                "| "
                f"[`{method.name}`](#{method_anchor})"
                f" | {table_cell(first_sentence(method.docstring))} |"
            )
            detail_blocks.append(render_python_method(symbol.name, method, method_anchor))
        method_index = "\n\n**Public Methods**\n\n" + "\n".join(rows)
        method_details = "\n\n" + "\n".join(detail_blocks)

    return f"""### `{symbol.name}` {{#{anchor}}}

{docstring_to_markdown(symbol.docstring)}
{method_index}
{method_details}
"""


def render_python_method(class_name: str, method: MethodDoc, anchor: str) -> str:
    return f"""#### `{class_name}.{method.name}` {{#{anchor}}}

<details className="api-detail">
<summary><code>{html.escape(class_name)}.{html.escape(method.name)}</code></summary>

{render_python_signature_block(method.signature)}

{docstring_to_markdown(method.docstring)}

</details>
"""


def render_python_signature_block(signature: str) -> str:
    return f"""<div className="api-signature">

```python
{signature}
```

</div>"""


def render_r_section(page: LoaderPage, r_docs: dict[str, RTopicDoc]) -> str:
    rows = ["| Topic | Summary |", "| --- | --- |"]
    details = []
    for topic_name in page.r_topics:
        topic = r_docs[topic_name]
        rows.append(
            "| "
            f"[`{topic_name}`](#{r_anchor_for_page(page, topic_name)})"
            f" | {table_cell(first_sentence(topic.description or topic.title))} |"
        )
        details.append(render_r_topic(topic, topic_name, r_anchor_for_page(page, topic_name)))
    return f"""## {page.title}

{page.summary}

{chr(10).join(rows)}

{chr(10).join(details)}
"""


def render_common_r_section(r_docs: dict[str, RTopicDoc]) -> str:
    rows = ["| Topic | Summary |", "| --- | --- |"]
    details = []
    for topic_name in COMMON_R_TOPICS:
        topic = r_docs[topic_name]
        anchor = common_r_anchor(topic_name)
        rows.append(
            "| "
            f"[`{topic_name}`](#{anchor})"
            f" | {table_cell(first_sentence(topic.description or topic.title))} |"
        )
        details.append(render_r_topic(topic, topic_name, anchor))
    return f"""## Common

Shared R loader API entries that do not belong to one output type.

{chr(10).join(rows)}

{chr(10).join(details)}
"""


def render_r_topic(topic: RTopicDoc, display_name: str, anchor: str) -> str:
    parts = [
        f"#### `{display_name}` {{#{anchor}}}",
        "",
        '<details className="api-detail">',
        f"<summary><code>{html.escape(display_name)}</code></summary>",
        "",
    ]
    if topic.description:
        parts.extend([escape_mdx_text(topic.description), ""])
    parts.extend(["```r", topic.usage.strip(), "```", ""])
    if topic.arguments:
        parts.extend(["**Arguments**", ""])
        for argument in topic.arguments:
            parts.append(f"- `{argument.name}`: {escape_mdx_text(argument.description)}")
        parts.append("")
    if topic.value:
        parts.extend(["**Returns**", "", escape_mdx_text(topic.value), ""])
    if topic.details:
        parts.extend(["**Details**", "", escape_mdx_text(topic.details), ""])
    if topic.examples:
        parts.extend(["**Examples**", "", "```r", topic.examples.strip(), "```", ""])
    parts.append("</details>")
    return "\n".join(parts)


def load_python_docs(
    package_dir: Path,
    public_symbols: set[str],
) -> dict[str, PythonSymbolDoc]:
    module_docs: dict[str, PythonSymbolDoc] = {}
    for module_path in sorted(package_dir.glob("*.py")):
        if module_path.name.startswith("_"):
            continue
        tree = ast.parse(module_path.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                module_docs[node.name] = PythonSymbolDoc(
                    name=node.name,
                    kind="function",
                    signature=render_python_signature(node, node.name),
                    docstring=ast.get_docstring(node) or "",
                    methods=(),
                )
            elif isinstance(node, ast.ClassDef):
                methods = tuple(
                    MethodDoc(
                        name=method.name,
                        signature=render_python_signature(
                            method,
                            f"{node.name}.{method.name}",
                            skip_first_argument=True,
                        ),
                        docstring=ast.get_docstring(method) or "",
                    )
                    for method in node.body
                    if isinstance(method, ast.FunctionDef) and is_public_method(method.name)
                )
                module_docs[node.name] = PythonSymbolDoc(
                    name=node.name,
                    kind="class",
                    signature=f"class {node.name}",
                    docstring=ast.get_docstring(node) or "",
                    methods=methods,
                )
    return {name: module_docs[name] for name in public_symbols if name in module_docs}


def parse_python_all(init_path: Path) -> tuple[str, ...]:
    tree = ast.parse(init_path.read_text(encoding="utf-8"))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets):
            continue
        values = ast.literal_eval(node.value)
        return tuple(str(value) for value in values if not str(value).startswith("__"))
    raise ValueError(f"Could not find __all__ in {init_path}")


def is_public_method(name: str) -> bool:
    return not name.startswith("_") and name not in {"__init__", "__repr__"}


def render_python_signature(
    node: ast.FunctionDef,
    display_name: str,
    *,
    skip_first_argument: bool = False,
) -> str:
    positional_arguments = list(node.args.posonlyargs) + list(node.args.args)
    if skip_first_argument and positional_arguments:
        positional_arguments = positional_arguments[1:]

    rendered_arguments: list[str] = []
    default_offset = len(positional_arguments) - len(node.args.defaults)
    for argument_index, argument in enumerate(positional_arguments):
        default_value = None
        if argument_index >= default_offset and node.args.defaults:
            default_value = node.args.defaults[argument_index - default_offset]
        rendered_arguments.append(render_python_argument(argument, default_value))

    if node.args.vararg:
        rendered_arguments.append("*" + render_python_argument(node.args.vararg, None))
    elif node.args.kwonlyargs:
        rendered_arguments.append("*")

    for argument, default_value in zip(node.args.kwonlyargs, node.args.kw_defaults):
        rendered_arguments.append(render_python_argument(argument, default_value))

    if node.args.kwarg:
        rendered_arguments.append("**" + render_python_argument(node.args.kwarg, None))

    return_annotation = ""
    if node.returns is not None:
        return_annotation = f" -> {ast.unparse(node.returns)}"

    return f"{display_name}({', '.join(rendered_arguments)}){return_annotation}"


def render_python_argument(argument: ast.arg, default_value: ast.expr | None) -> str:
    rendered = argument.arg
    if argument.annotation is not None:
        rendered += f": {ast.unparse(argument.annotation)}"
    if default_value is not None:
        rendered += f" = {ast.unparse(default_value)}"
    return rendered


def docstring_to_markdown(docstring: str) -> str:
    if not docstring.strip():
        return "_No docstring available._"
    lines = docstring.strip().splitlines()
    output: list[str] = []
    line_index = 0
    while line_index < len(lines):
        line = lines[line_index].rstrip()
        next_line = lines[line_index + 1].strip() if line_index + 1 < len(lines) else ""
        if is_numpy_docstring_section(lines, line_index):
            section_title = line.strip()
            line_index += 2
            section_lines: list[str] = []
            while line_index < len(lines) and not is_numpy_docstring_section(lines, line_index):
                section_lines.append(lines[line_index])
                line_index += 1
            output.append(render_numpy_docstring_section(section_title, section_lines))
            output.append("")
            continue
        output.append(line)
        line_index += 1
    return escape_mdx_text("\n".join(output).strip())


def is_numpy_docstring_section(lines: list[str], line_index: int) -> bool:
    section_titles = {"Parameters", "Returns", "Raises", "Examples"}
    if line_index + 1 >= len(lines):
        return False
    title = lines[line_index].strip()
    underline = lines[line_index + 1].strip()
    return title in section_titles and bool(underline) and set(underline) == {"-"}


def render_numpy_docstring_section(section_title: str, lines: list[str]) -> str:
    if section_title in {"Parameters", "Returns", "Raises"}:
        entries = parse_numpy_docstring_entries(lines)
        if not entries:
            return f"**{section_title}**"
        rendered_entries = []
        for entry_name, entry_description in entries:
            if entry_description:
                rendered_entries.append(f"- `{entry_name}`: {entry_description}")
            else:
                rendered_entries.append(f"- `{entry_name}`")
        return f"**{section_title}**\n\n" + "\n".join(rendered_entries)

    content = normalize_blank_lines("\n".join(lines))
    if not content:
        return f"**{section_title}**"
    return f"**{section_title}**\n\n```python\n{content}\n```"


def parse_numpy_docstring_entries(lines: list[str]) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    current_name = ""
    current_description_lines: list[str] = []
    for line in lines:
        if not line.strip():
            continue
        if line == line.lstrip():
            if current_name:
                entries.append(
                    (current_name, normalize_inline_whitespace(" ".join(current_description_lines)))
                )
            current_name = normalize_numpy_entry_name(line.strip())
            current_description_lines = []
            continue
        current_description_lines.append(line.strip())
    if current_name:
        entries.append(
            (current_name, normalize_inline_whitespace(" ".join(current_description_lines)))
        )
    return entries


def normalize_numpy_entry_name(name: str) -> str:
    if " : " not in name:
        return name
    entry_name, entry_type = name.split(" : ", 1)
    return f"{entry_name} ({entry_type})"


def load_r_docs(
    package_dir: Path,
    exported_names: set[str],
) -> dict[str, RTopicDoc]:
    topics: dict[str, RTopicDoc] = {}
    for rd_path in sorted((package_dir / "man").glob("*.Rd")):
        topic = parse_rd_topic(rd_path)
        if any(alias in exported_names for alias in topic.aliases):
            topics[topic.name] = topic
            for alias in topic.aliases:
                if alias in exported_names:
                    topics.setdefault(alias, topic)
    return topics


def parse_r_exports(namespace_path: Path) -> set[str]:
    exports = set()
    for line in namespace_path.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r"export\(([^)]+)\)", line.strip())
        if match:
            exports.add(match.group(1))
    return exports


def parse_rd_topic(rd_path: Path) -> RTopicDoc:
    text = "\n".join(
        line for line in rd_path.read_text(encoding="utf-8").splitlines() if not line.startswith("%")
    )
    name = rd_text(extract_first_rd_command(text, "name"))
    aliases = tuple(rd_text(alias) for alias in extract_rd_commands(text, "alias"))
    arguments = tuple(parse_rd_arguments(extract_first_rd_command(text, "arguments")))
    return RTopicDoc(
        name=name,
        aliases=aliases,
        title=rd_text(extract_first_rd_command(text, "title")),
        usage=rd_usage_text(extract_first_rd_command(text, "usage")),
        arguments=arguments,
        value=rd_text(extract_first_rd_command(text, "value")),
        description=rd_text(extract_first_rd_command(text, "description")),
        details=rd_text(extract_first_rd_command(text, "details")),
        examples=rd_examples_text(extract_first_rd_command(text, "examples")),
    )


def extract_first_rd_command(text: str, command_name: str) -> str:
    blocks = extract_rd_commands(text, command_name)
    return blocks[0] if blocks else ""


def extract_rd_commands(text: str, command_name: str) -> list[str]:
    blocks: list[str] = []
    search_start = 0
    command_prefix = f"\\{command_name}"
    while True:
        command_index = text.find(command_prefix, search_start)
        if command_index == -1:
            return blocks
        brace_index = text.find("{", command_index + len(command_prefix))
        if brace_index == -1:
            return blocks
        content, end_index = extract_braced_content(text, brace_index)
        blocks.append(content)
        search_start = end_index


def extract_braced_content(text: str, opening_brace_index: int) -> tuple[str, int]:
    depth = 0
    content_start = opening_brace_index + 1
    index = opening_brace_index
    while index < len(text):
        character = text[index]
        if character == "\\":
            index += 2
            continue
        if character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth == 0:
                return text[content_start:index], index + 1
        index += 1
    raise ValueError("Unclosed Rd brace block")


def parse_rd_arguments(arguments_block: str) -> list[RArgumentDoc]:
    arguments: list[RArgumentDoc] = []
    search_start = 0
    while True:
        item_index = arguments_block.find("\\item", search_start)
        if item_index == -1:
            return arguments
        name_brace_index = arguments_block.find("{", item_index)
        if name_brace_index == -1:
            return arguments
        name, name_end = extract_braced_content(arguments_block, name_brace_index)
        description_brace_index = arguments_block.find("{", name_end)
        if description_brace_index == -1:
            return arguments
        description, description_end = extract_braced_content(arguments_block, description_brace_index)
        arguments.append(
            RArgumentDoc(
                name=rd_text(name),
                description=rd_text(description),
            )
        )
        search_start = description_end


def rd_usage_text(text: str) -> str:
    return normalize_blank_lines(rd_text(text, preserve_newlines=True))


def rd_examples_text(text: str) -> str:
    if not text.strip():
        return ""
    unwrapped = unwrap_rd_examples(text)
    return normalize_blank_lines(rd_text(unwrapped, preserve_newlines=True))


def unwrap_rd_examples(text: str) -> str:
    text = text.strip()
    for command_name in ("dontrun", "donttest", "dontshow"):
        command = f"\\{command_name}"
        if text.startswith(command):
            brace_index = text.find("{", len(command))
            if brace_index != -1:
                content, _ = extract_braced_content(text, brace_index)
                return content
    return text


def rd_text(text: str, *, preserve_newlines: bool = False) -> str:
    if not text:
        return ""
    converted = text
    converted = convert_rd_itemize(converted)
    converted = re.sub(r"\\method\{([^{}]+)\}\{([^{}]+)\}", r"\1.\2", converted)
    for command_name, replacement in (
        ("code", r"`\1`"),
        ("verb", r"`\1`"),
        ("file", r"`\1`"),
        ("pkg", r"`\1`"),
        ("link", r"`\1`"),
        ("strong", r"**\1**"),
        ("emph", r"*\1*"),
    ):
        pattern = re.compile(rf"\\{command_name}\{{([^{{}}]*)\}}")
        previous = None
        while previous != converted:
            previous = converted
            converted = pattern.sub(replacement, converted)
    converted = re.sub(r"\\url\{([^{}]+)\}", r"\1", converted)
    converted = converted.replace("\\cr", "\n")
    converted = converted.replace("\\%", "%")
    converted = converted.replace("\\_", "_")
    converted = converted.replace("\\{", "{")
    converted = converted.replace("\\}", "}")
    converted = converted.replace("\\", "")
    if preserve_newlines:
        return "\n".join(line.rstrip() for line in converted.strip().splitlines())
    return normalize_inline_whitespace(converted)


def convert_rd_itemize(text: str) -> str:
    converted = text
    search_start = 0
    while True:
        itemize_index = converted.find("\\itemize", search_start)
        if itemize_index == -1:
            return converted
        brace_index = converted.find("{", itemize_index + len("\\itemize"))
        if brace_index == -1:
            return converted
        content, end_index = extract_braced_content(converted, brace_index)
        raw_items = [item.strip() for item in re.split(r"\\item\s+", content) if item.strip()]
        replacement = " ".join(raw_items)
        converted = converted[:itemize_index] + replacement + converted[end_index:]
        search_start = itemize_index + len(replacement)


def normalize_inline_whitespace(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def normalize_blank_lines(text: str) -> str:
    lines = [line.rstrip() for line in text.strip().splitlines()]
    output: list[str] = []
    previous_blank = False
    for line in lines:
        is_blank = not line.strip()
        if is_blank and previous_blank:
            continue
        output.append(line)
        previous_blank = is_blank
    return "\n".join(output).strip()


def validate_loader_config(
    python_public_symbols: set[str],
    python_docs: dict[str, PythonSymbolDoc],
    r_public_topics: set[str],
    r_docs: dict[str, RTopicDoc],
) -> None:
    configured_python_symbols = assigned_python_symbols()
    configured_r_topics = assigned_r_topics()

    messages = []
    messages.extend(
        validate_configured_symbols_exist(
            "configured Python symbols without source docs",
            configured_python_symbols,
            set(python_docs),
            "LOADER_PAGES or COMMON_PYTHON_SYMBOLS",
        )
    )
    messages.extend(
        validate_configured_symbols_exist(
            "configured R topics without roxygen docs",
            configured_r_topics,
            set(r_docs),
            "LOADER_PAGES or COMMON_R_TOPICS",
        )
    )
    messages.extend(
        validate_public_symbols_are_assigned(
            "public Python symbols not assigned to loader docs",
            python_public_symbols,
            configured_python_symbols,
            "LOADER_PAGES or COMMON_PYTHON_SYMBOLS",
        )
    )
    messages.extend(
        validate_public_symbols_are_assigned(
            "public R exports not assigned to loader docs",
            r_public_topics,
            configured_r_topics,
            "LOADER_PAGES or COMMON_R_TOPICS",
        )
    )
    messages.extend(validate_python_docstrings(python_public_symbols, python_docs))
    messages.extend(validate_r_docstrings(r_public_topics, r_docs))

    if messages:
        joined_messages = "\n".join(f"- {message}" for message in messages)
        raise ValueError(f"Loader docs configuration drift:\n{joined_messages}")


def assigned_python_symbols() -> set[str]:
    symbols = set(COMMON_PYTHON_SYMBOLS)
    for page in LOADER_PAGES:
        symbols.update(page.python_symbols)
    return symbols


def assigned_r_topics() -> set[str]:
    topics = set(COMMON_R_TOPICS)
    for page in LOADER_PAGES:
        topics.update(page.r_topics)
    return topics


def validate_configured_symbols_exist(
    label: str,
    configured_symbols: set[str],
    documented_symbols: set[str],
    config_location: str,
) -> list[str]:
    missing_symbols = sorted(configured_symbols - documented_symbols)
    if not missing_symbols:
        return []
    return [
        f"{label}: {', '.join(missing_symbols)}. "
        f"Remove them from {config_location} or add source documentation."
    ]


def validate_public_symbols_are_assigned(
    label: str,
    public_symbols: set[str],
    configured_symbols: set[str],
    config_location: str,
) -> list[str]:
    unassigned_symbols = sorted(public_symbols - configured_symbols)
    if not unassigned_symbols:
        return []
    return [
        f"{label}: {', '.join(unassigned_symbols)}. "
        f"Assign them in {config_location} so the generated website stays current."
    ]


def validate_python_docstrings(
    python_public_symbols: set[str],
    python_docs: dict[str, PythonSymbolDoc],
) -> list[str]:
    missing_symbols = sorted(
        symbol_name
        for symbol_name in python_public_symbols
        if symbol_name not in python_docs or not python_docs[symbol_name].docstring.strip()
    )
    missing_methods = sorted(
        f"{symbol.name}.{method.name}"
        for symbol_name in python_public_symbols
        if symbol_name in python_docs
        for symbol in [python_docs[symbol_name]]
        for method in symbol.methods
        if not method.docstring.strip()
    )
    messages = []
    if missing_symbols:
        messages.append(
            "public Python symbols without docstrings: " + ", ".join(missing_symbols)
        )
    if missing_methods:
        messages.append(
            "public Python methods without docstrings: " + ", ".join(missing_methods)
        )
    return messages


def validate_r_docstrings(
    r_public_topics: set[str],
    r_docs: dict[str, RTopicDoc],
) -> list[str]:
    missing_topics = sorted(
        topic_name
        for topic_name in r_public_topics
        if topic_name not in r_docs
        or not (r_docs[topic_name].description or r_docs[topic_name].title).strip()
    )
    if not missing_topics:
        return []
    return ["public R exports without roxygen documentation: " + ", ".join(missing_topics)]


def command_doc_slug(command: str) -> str:
    return command.removeprefix("cfdna ").strip()


def first_sentence(text: str) -> str:
    text = normalize_inline_whitespace(text)
    if not text:
        return "Reference entry."
    match = re.search(r"(?<=[.?!])\s+", text)
    return text[: match.start()].strip() if match else text


def table_cell(text: str) -> str:
    return escape_mdx_text(text).replace("|", "\\|").replace("\n", " ")


def escape_mdx_text(text: str) -> str:
    return text.replace("{", "&#123;").replace("}", "&#125;")


def section_anchor(title: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")


def python_anchor(symbol_name: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_-]+", "-", symbol_name).lower()


def r_anchor(topic_name: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_-]+", "-", topic_name).lower()


def r_anchor_for_page(page: LoaderPage, topic_name: str) -> str:
    return f"{section_anchor(page.title)}-{r_anchor(topic_name)}"


def common_r_anchor(topic_name: str) -> str:
    return f"common-{r_anchor(topic_name)}"


if __name__ == "__main__":
    main()
