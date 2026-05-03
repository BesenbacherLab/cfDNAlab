# `cfdna bam-to-bam`

Apply cfDNAlab fragment filters and optional correction weights to an existing BAM file. The command writes the surviving original BAM records with fragment-level metadata attached as AUX tags.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    setup["Prepare run<br/>validate options, resolve chromosomes, load side files"]
    output_writer["Open output BAM<br/>with the input header"]

    input --> setup --> output_writer

    subgraph side_inputs["Optional side inputs"]
        windows["BED windows<br/>keep overlapping fragments"]
        blacklist["Blacklist BED<br/>exclude problematic regions"]
        scaling["Scaling TSVs<br/>coverage and fragment-count weights"]
        gc["GC package + 2bit reference<br/>sample-specific GC correction"]
    end

    windows --> setup
    blacklist --> setup
    scaling --> setup
    gc --> setup

    subgraph chromosome_pass["Chromosome pass"]
        fetch["Fetch selected reads<br/>whole chromosome or BED-focused span"]
        read_filter{"Read-level filters"}
        assemble["Build fragments<br/>paired-end or reads-as-fragments"]
        fragment_filter{"Fragment-level filters"}
        weight["Compute fragment metadata<br/>length, GC weight, scaling weights"]
        tag_records["Attach metadata<br/>to each surviving BAM record"]
        sort_records["Sort records within chromosome"]

        fetch --> read_filter
        read_filter -- "pass" --> assemble
        read_filter -- "fail" --> dropped_read["Drop read"]
        assemble --> fragment_filter
        fragment_filter -- "pass" --> weight --> tag_records --> sort_records
        fragment_filter -- "fail" --> dropped_fragment["Drop fragment"]
    end

    output_writer --> fetch
    sort_records --> final_bam["Filtered and tagged BAM"]
    final_bam --> stats["Run statistics<br/>fragments included, skipped, and filtered"]

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef decision fill:#fff5d9,stroke:#b88700,color:#2f2500;
    classDef drop fill:#fbe7e7,stroke:#b84545,color:#3a1111;
    classDef output fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,setup,output_writer,fetch,assemble,weight,tag_records,sort_records core;
    class windows,blacklist,scaling,gc optional;
    class read_filter,fragment_filter decision;
    class dropped_read,dropped_fragment drop;
    class final_bam,stats output;
```

## What Is Preserved

For each surviving fragment, `bam-to-bam` writes the original BAM record or records. It preserves flags, CIGARs, sequences, qualities, mate fields, and read names.

## What Is Added

Each surviving record receives the fragment length. When requested, records also receive GC correction weights, coverage-scaling weights, and fragment-count-scaling weights. Paired-end fragments write the same fragment-level metadata to both mates. Unpaired `--reads-are-fragments` mode writes one record per fragment.
