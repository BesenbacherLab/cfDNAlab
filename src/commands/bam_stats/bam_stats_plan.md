# BAM Stats Plan

## Purpose

- Inventory what is in a BAM file with clear, auditable summaries
- Provide fast QC and integrity checks for cfDNA pipelines
- Report both read-level and fragment-level summaries without feature extraction

## Non-goals and relationship to existing commands

- Do not recreate outputs from `cfdna lengths`, `cfdna fcoverage`, `cfdna fragment-kmers`,
  `cfdna transitions`, `cfdna wps`, `cfdna wps-peaks`, or `cfdna midpoints`
- Avoid windowed or per-base outputs by default
- Avoid reference or annotation requirements unless explicitly enabled
- Keep metrics small, interpretable, and easy to validate

## Output formats

- CLI summary for humans
- JSON for pipelines
- TSV for spreadsheets and quick plotting
- Optional per-read and per-fragment summary tables with stable column definitions

## Phases

### Phase 0: Define scope and outputs

- Agree on the first release set and what is optional or heavy
- Define the minimum report structure for every run
- Decide default histogram bins and thresholds

### Phase 1: Enumerate metrics and aggregation axes

- Build a matrix of metrics with required tags and prerequisites
- Specify if each metric is read-level, fragment-level, or both
- Identify edge cases and invalid input handling

### Phase 2: Data flow and performance plan

- Decide on single pass vs multi pass trade-offs
- Define memory ceilings and streaming aggregation
- Identify metrics that require mate lookups or full read pairing

### Phase 3: Validation and error model

- Decide strict vs permissive parsing
- Define warnings for malformed flags, missing mates, or invalid CIGAR
- Provide clear error messages for unsupported inputs

### Phase 4: Test plan

- Design tiny BAM fixtures with known outputs
- Create deterministic expectation tables per metric
- Add edge case fixtures with missing tags and invalid records

## Core aggregation axes

These axes apply to many metrics and reduce repetition.

- Overall
- By contig
- By read group
- By library and sample when present
- By MAPQ bins
- By read length bins
- By insert size bins
- By strand
- By read number for pairs
- By filter tier

## Filter strata and attrition

Define shared filter tiers and report every core metric for each tier. Also report a drop
table that shows how many records or fragments were removed by each step, with the most
common reasons. These tiers should reuse the same flag logic used across commands.

- Raw: no filtering beyond valid record decoding
- Primary only: exclude secondary and supplementary
- Mapped pair: read and mate mapped on same contig
- Orientation and pairing: inward orientation and proper pair if required
- QC pass: remove duplicate and QC failed records
- MAPQ pass: mapq >= configured threshold
- Fragment length pass: inside configured length bounds
- Default filter set: the shared include and fragment filters used by cfDNA commands

## Statistics catalog

### Core report

#### File metadata and integrity

- Header summary: contigs and lengths, sort order, read groups, program records
- Presence of index file and reported per-contig read counts from the index
- Compression format and version if available
- Records with missing SEQ or QUAL
- Records that reference contigs not present in the header

#### Record classification counts

- Total records
- Mapped, unmapped, primary, secondary, supplementary
- Paired, unpaired, mixed
- QC failed records
- Duplicate records
- Read1 vs read2 counts

#### Flag occurrences

Report counts and rates for every SAM flag at read and fragment levels.

- Read level: direct flag counts from each record
- Fragment level: derived counts from the pair
- Examples of fragment derivations:
  - Either read unmapped
  - Both reads mapped
  - Either read duplicate
  - Both reads QC failed
  - Proper pair when both reads present and aligned

#### Read name and pairing integrity

- Unique read name count
- Records per read name distribution: 1, 2, more than 2
- Read names with missing mate record
- Read names with multiple primary alignments
- Read names with conflicting read1 and read2 flags

#### Mapping and alignment summaries

- MAPQ distribution: min, max, mean, median, histogram
- MQ0 rate and fraction of reads with MAPQ below threshold
- Reference usage: counts and fractions by contig
- Read orientation by contig and overall
- Alignment length distribution from CIGAR
- Reference span distribution
- Fraction of reads with spliced alignments based on N in CIGAR

#### CIGAR operations and clipping

- Totals and rates of each CIGAR op: M, I, D, N, S, H, P, =, X
- Soft clipping and hard clipping counts and sizes
- Clipping location: 5' vs 3' vs both
- Clipping rate by contig and by read group
- Insertions and deletions: counts, lengths, rate per aligned base
- CIGAR operation count per read as a proxy for alignment complexity

#### Read sequence and quality summaries

- Read length distribution from SEQ and from CIGAR
- Per-read mean, min, max base quality
- Aggregate base quality histogram
- Low quality base counts based on threshold
- Base composition and GC content per read
- Fraction of reads with ambiguous bases

#### Tag presence and basic validity

- Tag presence counts for all observed tags
- RG tag coverage and unknown RG counts
- UMI and barcode tag presence when present
- Presence of NM, MD, AS, SA, MC, and MQ tags

#### Fragment and pairing summaries

- Proper pair fraction
- Mate orientation: FR, RF, FF, RR
- Insert size summary from TLEN: min, max, mean, median, quantiles, coarse histogram
- Outlier fraction for insert size by thresholds
- Interchromosomal pairs and same-contig pairs
- Fragment length summary from positions when TLEN is unreliable
- Overlap length and overlap fraction for paired reads
- Gap size summary for non-overlapping pairs
- TLEN sign distribution and discordant sign rates
- Read1 and read2 strand balance
- Inward vs outward oriented pairs
- Fraction of pairs that pass the cfDNA fragment definition used by tools

#### cfDNA fragmentomics overview

These are coarse summaries for inventory, not full feature extraction.

- Fragment length quantiles and coarse bins
- Mono, di, and tri nucleosome window counts and ratios
- Short fragment rate below a minimum threshold
- Long fragment tail rate above a maximum threshold
- Fragment length periodicity index from the length histogram
- Canonical cfDNA band fraction for a user-defined band

#### Error and sanity checks

- Missing mate for paired reads
- Inconsistent flags such as proper pair with unmapped
- Invalid or empty CIGAR
- Negative or zero-length reference spans
- Mates mapped to unknown contigs
- TLEN outliers when both reads are mapped to the same contig
- Records that violate declared sort order
- RG tags not declared in the header
- Paired reads on different contigs with nonzero TLEN

### Expanded report

#### Extended tag summaries

- Tag type mismatch counts, for example integer vs string for the same tag
- NM tag edit distance summaries
- MD tag mismatch summaries when present
- AS tag alignment score summaries
- SA tag count per read name
- MC tag presence and mate CIGAR length summaries
- Tag value range checks for known tags

#### Multi-alignment and split-read summaries

- Secondary count per read name
- Supplementary count per read name
- Fraction of reads with multiple alignments
- Split read patterns: number of segments per read

#### Expanded cfDNA fragmentomics

- Fragment length by read group and by MAPQ bins
- Fragment length by contig and by strand
- Fragment length by coarse GC bins
- Periodicity index by read group
- Mono to di nucleosome ratio by read group
- Short and long fragment rates by read group

#### Optional heavy summaries

Keep disabled by default and make cost visible.

- Per-cycle base composition and quality
- Read start and read end base composition summaries
- Coverage depth summaries by contig, without per-base outputs

## Notes on scope

- The command is for an inventory of the BAM file, not a replacement for existing feature tools
- Keep results stable and easy to compare across runs
- Add new metrics only if they improve understanding of what is inside the BAM
