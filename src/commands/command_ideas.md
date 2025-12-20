Generated ideas for future commands.

## Fragmentomics-like feature ideas not covered today
- Fragment end coverage for 5' ends only, with strand-aware counts
- End site diversity per window, using entropy or Gini to describe end clustering
- End site spacing distribution and periodicity from adjacent end positions
- Fragment end orientation asymmetry around genomic landmarks
- Fragment end enrichment at open chromatin or promoter regions, when annotations are provided
- Fragment end hotspot detection and enrichment ratios vs local background
- Fragment end position jaggedness, including 1 bp and 10 bp offset patterns
- Fragment end distance to feature boundaries, such as TSS or TFBS
- Fragment end base composition by position around the cut site, using reference bases
- Fragment end phasing scores from periodicity in cut site density
- Fragment end skew around gene bodies, stratified by gene strand
- Preferred end sites, defined as high frequency cut sites within windows
- End site dispersion within windows, using median absolute deviation
- Fragment end distance to nearest nucleosome peak, when peaks are provided
- End site density around copy number breakpoints, when intervals are provided
- Paired end motif co-occurrence, stratified by fragment length
- End motif asymmetry by strand
- End position enrichment around nucleosome peaks from `wps-peaks` outputs
- End phase shift relative to WPS or nucleosome peak centers

## Derived features from existing outputs
These can be computed from the windowed output of `cfdna lengths` and should not be new commands.
- Short to long fragment ratios inside coarse genomic bins
- Size band ratios for canonical cfDNA bands, reported by contig or read group
- Fragment length quantiles, mean, and variance by window
- Fragment length periodicity index from the length histogram
- Mono, di, and tri nucleosome band ratios from length bins

## Less obvious tumor or tissue signal ideas
These are higher risk and higher reward ideas that could help classify cancer vs control
or estimate tissue of origin. Most need optional annotations or external atlases.
- Fragment end enrichment against tissue specific open chromatin atlases
- End density phase shift relative to known nucleosome dyads in different tissues
- Cut site periodicity strength differences between active and inactive chromatin
- End motif ratios tied to nuclease activity signatures, such as DNASE1L3 vs DFFB
- Motif to length coupling, where motif frequencies shift by fragment size bands
- Strand aware end imbalance around TSS, TES, and enhancers by gene strand
- End clustering score near TF binding sites using motif catalogs or ATAC peaks
- Short fragment enrichment in promoter regions vs gene bodies
- Fragment length skew within CpG islands vs flanking regions
- End site entropy differences between active and repressed chromatin states
- End hotspot concordance with known cancer breakpoints or CNV segments
- Fragment end density changes near replication timing domains
- End site distance to nucleosome peaks from `wps-peaks` and external maps
- Fragment size variance near histone mark peaks, such as H3K27ac
- End density imbalance across compartments, such as A and B compartments
- Mitochondrial read fraction and fragment length skew in mtDNA
- End site asymmetry around immune cell specific enhancers for tissue deconvolution
- Fragmentation signatures around tissue specific methylation proxy regions

## Command ideas to make these available
- `cfdna end-hotspots`:
  - Identify windows with enriched end density over local background
  - Output scores and ranked windows for downstream review
- `cfdna end-entropy`:
  - Compute entropy or Gini per window from end counts
  - Surface highly clustered or overly uniform fragmentation patterns
- `cfdna end-spacing`:
  - Compute spacing distributions between adjacent ends
  - Report periodicity indices and peak spacing summaries
- `cfdna end-orientation`:
  - Summarize forward and reverse end imbalance by window or feature
  - Report asymmetry ratios and their distributions
- `cfdna end-context`:
  - Summarize reference base composition around cut sites
  - Report position specific base frequencies and GC composition
- `cfdna end-periodicity`:
  - Compute periodicity indices from cut site density
  - Report dominant periods and strength scores
- `cfdna end-dispersion`:
  - Compute dispersion of cut sites per window
  - Flag windows with unusually clustered ends
- `cfdna end-to-feature`:
  - Summarize distances from cut sites to annotated features
  - Report distributions by feature class
- `cfdna end-pairs`:
  - Count paired end motif combinations for both fragment ends
  - Stratify by fragment length bins
- `cfdna end-to-peaks`:
  - Summarize distances from cut sites to nucleosome peaks
  - Accept peaks from `wps-peaks` or a BED file
- `cfdna end-footprints`:
  - Compute summary scores around annotated sites without full profiles
  - Report dip depth, flank height, and phasing amplitude
