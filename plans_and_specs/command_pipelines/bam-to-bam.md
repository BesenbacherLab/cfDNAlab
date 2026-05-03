# `cfdna bam-to-bam` pipeline

Purpose: read an indexed coordinate-sorted BAM, apply cfDNAlab fragment filters and optional correction weights, and write a BAM containing the surviving original records with fragment-level AUX metadata.

## Flow

```mermaid
flowchart TD
    input["Input BAM<br/>plus BAI/CSI index"]
    validate["Validate command configuration<br/>fragment lengths, GC reference, paired/unpaired constraints"]
    resolve["Resolve selected chromosomes<br/>and BAM contig lengths"]
    load["Load optional run-wide inputs<br/>BED windows, blacklists, scaling TSVs, GC package, 2bit reference"]
    writer["Open output BAM writer<br/>using the input BAM header"]

    input --> validate --> resolve --> load --> writer

    subgraph per_chromosome["Per selected chromosome"]
        chrom_reader["Open indexed chromosome reader"]
        fetch["Fetch reads<br/>whole chromosome or BED span plus fragment-length halo"]
        stream["Stream records<br/>in coordinate order"]
        read_filter{"Read passes<br/>always-on filters?"}
        skip_read["Drop read"]
        fragment_mode{"Fragment mode"}
        paired["Paired-end mode<br/>pair by qname, require inward orientation<br/>span = forward.pos to reverse.reference_end"]
        unpaired["Unpaired mode<br/>one read is one fragment<br/>span = read.pos to read.reference_end"]
        fragment_filter{"Fragment passes filters?<br/>length, blacklist, optional BED overlap"}
        skip_fragment["Drop fragment"]
        weights["Compute optional fragment weights<br/>GC, coverage scaling, count scaling"]
        tags["Attach fragment-level tag payload<br/>to each surviving BAM record"]
        sort_buffer["Bounded coordinate sorter<br/>orders records within the chromosome"]
        flush["Flush chromosome tail"]

        chrom_reader --> fetch --> stream --> read_filter
        read_filter -- "no" --> skip_read
        read_filter -- "yes" --> fragment_mode
        fragment_mode -- "paired BAM" --> paired
        fragment_mode -- "reads are fragments" --> unpaired
        paired --> fragment_filter
        unpaired --> fragment_filter
        fragment_filter -- "no" --> skip_fragment
        fragment_filter -- "yes" --> weights --> tags --> sort_buffer --> flush
    end

    writer --> chrom_reader
    flush --> output["Output BAM<br/>plus terminal/log fragment statistics"]

    sorted_warning["Review note B2B-001<br/>chromosome processing order can break BAM header target-id sortedness"]
    writer -.-> sorted_warning

    tag_warning["Review note G-017<br/>public AUX tag names need two-byte BAM-safe names"]
    tags -.-> tag_warning

    classDef inputOutput fill:#edf6ff,stroke:#3b73b9,color:#10233f;
    classDef process fill:#f7f7f4,stroke:#777,color:#202020;
    classDef decision fill:#fff4d6,stroke:#b88700,color:#2f2500;
    classDef drop fill:#fbe7e7,stroke:#b84545,color:#3a1111;
    classDef warning fill:#f3e8ff,stroke:#7c3fb2,color:#241033;

    class input,output inputOutput;
    class validate,resolve,load,writer,chrom_reader,fetch,stream,paired,unpaired,weights,tags,sort_buffer,flush process;
    class read_filter,fragment_mode,fragment_filter decision;
    class skip_read,skip_fragment drop;
    class sorted_warning,tag_warning warning;
```

## Notes

`bam-to-bam` keeps the original BAM records, flags, CIGARs, sequences, qualities, and mate fields for surviving fragments. Paired-end fragments write both mates; unpaired `--reads-are-fragments` mode writes the single read.

The command currently preserves the input BAM header while choosing its own chromosome processing order. Review finding B2B-001 tracks the resulting coordinate-sort risk.

The documented multi-character AUX tag names are under review in G-017. BAM AUX tag keys are two bytes, so the public tag vocabulary needs to be made explicit before release.
