# `cfdna bam-to-frag`

Convert a BAM file into a compressed fragment table. Each output row represents one fragment, with optional correction weights written as extra columns and described by a companion header file.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    reads["Read alignments<br/>from selected chromosomes"]
    read_filter{"Read-level filters"}
    assemble["Create fragments<br/>pair inward-facing mates or use each read as a fragment"]
    fragment_filter{"Fragment filters"}
    metadata["Compute row metadata<br/>min MAPQ, read1 strand, optional weights"]
    sort_rows["Sort fragment rows<br/>within each chromosome"]
    chunks["Write chromosome chunks"]
    merge["Merge chunks<br/>in chromosome order"]
    frag["Compressed fragment table<br/>frag.tsv.gz"]
    header["Column header file<br/>frag.header.tsv"]
    stats["Run statistics<br/>fragments included, skipped, and filtered"]

    filter_inputs["Optional filter inputs<br/>BED windows and blacklists"]
    weight_inputs["Optional weight inputs<br/>GC package, 2bit reference, scaling TSVs"]

    input --> reads --> read_filter
    read_filter -- "pass" --> assemble --> fragment_filter
    read_filter -- "fail" --> dropped_read["Drop read"]
    fragment_filter -- "pass" --> metadata --> sort_rows --> chunks --> merge --> frag --> stats
    fragment_filter -- "fail" --> dropped_fragment["Drop fragment"]
    metadata --> header

    filter_inputs -.-> fragment_filter
    weight_inputs -.-> metadata

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef decision fill:#fff5d9,stroke:#b88700,color:#2f2500;
    classDef drop fill:#fbe7e7,stroke:#b84545,color:#3a1111;
    classDef output fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,reads,assemble,metadata,sort_rows,chunks,merge core;
    class filter_inputs,weight_inputs optional;
    class read_filter,fragment_filter decision;
    class dropped_read,dropped_fragment drop;
    class frag,header,stats output;
```

## Fragment Rows

The base columns are chromosome, start, end, minimum MAPQ, and read1 strand. For paired-end input, start and end come from the inward-facing fragment span. In `--reads-are-fragments` mode, each accepted read becomes one fragment row.

## Optional Metadata

When requested, `bam-to-frag` adds GC correction weights, coverage-scaling weights, and fragment-count-scaling weights as extra columns. The companion `frag.header.tsv` file records the exact column names so downstream commands can restore the metadata.
