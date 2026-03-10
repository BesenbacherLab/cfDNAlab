<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna frag-to-bam

Convert a finaleDB-style frag file to a BAM file with unpaired reads (each read is a full fragment).

Each read in the new BAM file represents a fragment from the frag file.

The first five columns in the frag file should be: `Chromosome, Start, End, MapQ, Strand`.

## Extra columns

Optional extra columns can be transferred to BAM AUX tags when column names are known.

The recognized names and respective AUX tags are:

- `gc_weight` -> `GC`

- `scaling_weight` -> `COV`

- `flen` -> `FLEN`

## BAM file

The BAM header contains all contigs from `--chrom-sizes` in the `--chrom-sizes` order.

The BAM file is not indexed. This can be done with `samtools index`.


<hr class="cli-usage-separator" />

## Usage

`cfdna frag-to-bam [OPTIONS] --frag <FRAG> --output-dir <OUTPUT_DIR> --chrom-sizes <CHROM_SIZES>`

## Options

- `-h, --help`

  Print help (see a summary with '-h')
  

## Core

- `-i, --frag <FRAG>`

  Path to a coordinate-sorted `.tsv` frag file `[path]`
  

- `-o, --output-dir <OUTPUT_DIR>`

  Output directory to write new BAM file in `[path]`
  

- `-x, --output-prefix <OUTPUT_PREFIX>`

  Prefix for output file (e.g., a sample name) `[string]`
  
  E.g., specify to enable writing to the same output directory from multiple calls to this software.
  
  Examples produce files like: `<prefix>.bam`,
  
  [default: fragments]
  

- `--frag-header <FRAG_HEADER>`

  Optional header file with tab-separated column names for the frag file [path]
  
  Supply this when you want to transfer extra columns (`gc_weight`, `scaling_weight`, and/or `flen`) to AUX tags in the BAM file and the frag file has no inline header row.
  
  **Auto-detection**: The command also tries to auto-detect a companion header file named `<prefix>.frag.header.tsv` when the frag path follows `<prefix>.frag.tsv` (optionally with `.gz` or `.zst`).
  
  When no headers are supplied/detected or found inline, the command still accepts headerless 5-column frag files.
  
  Use `--ignore-extras` when you want to ignore all extra columns after the first five.
  

- `--chrom-sizes <CHROM_SIZES>`

  File with chromosome sizes (FAI or two-column sizes) for the BAM header `[path]`
  
  E.g. the UCSC `hg38.chrom.sizes` file (or similar for your assembly).
  

## Chromosome Selection (select max. one arg.)

- `--chromosomes <CHROMOSOMES>...`

  Names of chromosomes to process (comma-separated or repeated). E.g. `'chr1,chr2,chr3'`.
  
  When no chromosomes are specified, it defaults to `chr1..chr22`.
  
  Specify `"all"` *as the only string* to use all present chromosomes. For BAM-backed commands this uses the BAM header order. For commands that read chromosome order from their input, this may use the input order or some other order.
  

- `--chromosomes-file <CHROMOSOMES_FILE>`

  File with chromosome names to process (one per line)
  

## Filtering

- `--min-fragment-length <MIN_FRAGMENT_LENGTH>`

  Minimum fragment length to include `[integer]`
  
  [default: 30]
  

- `--max-fragment-length <MAX_FRAGMENT_LENGTH>`

  Maximum fragment length to include `[integer]`
  
  [default: 1000]
  

- `--min-mapq <MIN_MAPQ>`

  Minimum mapping quality to include `[integer]`
  
  Defaults to 0 to allow making filtering decisions downstream.
  
  [default: 0]
  

- `--ignore-extras`

  Ignore all frag columns after the first five `[flag]`
  
  This disables mapping extra columns to BAM AUX tags. It also allows headers with extra names that are not supported for AUX mapping.
  

- `--allow-unknown-extras`

  Allow unknown extra header columns and ignore them `[flag]`
  
  By default, unknown extra columns cause an error to prevent silent mistakes.
  
  With this flag, unknown extra columns are ignored with a warning, while known extra columns (`gc_weight`, `scaling_weight`, `flen`) are still transferred.
  
  If you want to ignore all extras, use `--ignore-extras` instead.
  

- `-b, --blacklist <BLACKLIST>...`

  Optional BED file(s) with blacklisted regions `[path]`
  

- `--blacklist-min-size <BLACKLIST_MIN_SIZE>`

  Minimum size of blacklist intervals to load (bp) `[integer]`
  
  [default: 1]
  

- `--blacklist-strategy <BLACKLIST_STRATEGY>`

  The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
  
  Possible values: `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
  
  Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
  
  [default: any]

