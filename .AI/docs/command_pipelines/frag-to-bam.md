# `cfdna frag-to-bam`

Convert a fragment table back into a BAM file. Each accepted frag row becomes one unpaired BAM record spanning the full fragment interval.

## Pipeline

```mermaid
flowchart TD
    inputs["Inputs<br/>fragment table, chrom sizes, optional column header"]
    layout["Resolve columns<br/>core fields and optional metadata"]
    parse["Read fragment rows<br/>chromosome, start, end, MAPQ, strand"]
    fragment_filter{"Fragment filters"}
    stage["Stage rows<br/>by chromosome"]
    header["Build BAM header<br/>from chrom sizes"]
    records["Create BAM records<br/>one unpaired record per fragment"]
    tags["Attach optional AUX metadata<br/>GC, scaling, fragment length"]
    write["Write BAM<br/>in chrom-sizes order"]
    bam["Fragment BAM<br/>fragments.bam"]
    stats["Run statistics<br/>parsed, filtered, and written fragments"]

    filter_inputs["Optional filter inputs<br/>chromosome selection and blacklists"]

    inputs --> layout --> parse --> fragment_filter
    fragment_filter -- "pass" --> stage --> header --> records --> tags --> write --> bam --> stats
    fragment_filter -- "fail" --> dropped_fragment["Drop fragment"]

    filter_inputs -.-> fragment_filter

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef decision fill:#fff5d9,stroke:#b88700,color:#2f2500;
    classDef drop fill:#fbe7e7,stroke:#b84545,color:#3a1111;
    classDef output fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class inputs,layout,parse,stage,header,records,tags,write core;
    class filter_inputs optional;
    class fragment_filter decision;
    class dropped_fragment drop;
    class bam,stats output;
```

## BAM Records

The output BAM contains one unpaired record per surviving fragment. The record starts at the frag `start`, ends at the frag `end`, uses a full-length match CIGAR, and stores the frag MAPQ and strand.

## Optional Metadata

If column names are available from an inline header, explicit header file, or companion `frag.header.tsv`, recognized extra columns are transferred to BAM AUX metadata. Headerless five-column files are also accepted when no extra metadata needs to be restored.
