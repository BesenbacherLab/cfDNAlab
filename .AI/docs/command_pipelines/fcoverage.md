# `cfdna fcoverage`

Compute fragment coverage from a BAM file. The command turns alignments into fragments, counts their covered bases, and writes either positional coverage or window-level summaries.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    tiles["Scan genome in tiles<br/>with fragment length halos"]
    fragments["Create fragments<br/>pair mates or use each read as a fragment"]
    filters["Keep usable fragments<br/>MAPQ, pairing, and length"]
    weights["Apply requested weights<br/>GC, scaling, and length normalization"]
    coverage["Accumulate coverage<br/>over counted fragment bases"]
    regions["Apply output regions<br/>windows, groups, and blacklists"]
    summarize["Choose coverage view<br/>per base, per window, or per group"]
    merge["Merge tile results<br/>and reduce windows crossing tile edges"]
    output["Coverage output<br/>bedGraph or compressed TSV"]
    group_index["Optional group index<br/>maps group_idx to group name"]
    stats["Run statistics<br/>fragments counted, skipped, and filtered"]

    region_inputs["Optional region inputs<br/>BED windows, grouped BED, fixed bins, blacklist"]
    correction_inputs["Optional correction inputs<br/>GC package, 2bit reference, scaling TSVs"]

    input --> tiles --> fragments --> filters --> weights --> coverage --> regions --> summarize --> merge --> output --> stats
    summarize --> group_index

    region_inputs -.-> tiles
    region_inputs -.-> regions
    correction_inputs -.-> weights

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,tiles,fragments,filters,weights,coverage,regions,summarize,merge core;
    class region_inputs,correction_inputs,group_index optional;
    class output,stats outputClass;
```

## Fragment Coverage Model

For paired-end BAM input, `fcoverage` builds one fragment from inward-facing mates and counts coverage across the fragment span. In `--reads-are-fragments` mode, each accepted read is counted as its own fragment. When gap-aware counting is enabled, deleted or skipped reference regions can be excluded from the counted bases.

## Output Views

Without windowing, the command writes per-position bedGraph coverage. With `--by-size` or `--by-bed`, it writes one compressed TSV row per window. With `--by-grouped-bed`, windows are reduced to group-level rows and a group-index TSV records the group names. BED mode can also write positional coverage restricted to the selected windows.

## Corrections And Weights

Optional GC correction, genome scaling, and length normalization modify each fragment's contribution before coverage is accumulated. Blacklists remove masked bases from window denominators and summary statistics while preserving their counts in aggregate outputs.
