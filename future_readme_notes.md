Info for readme for upcoming commands:




| `cfdna fragment-kmers`   | Count fragment k-mers with highly flexible selection of positions                                                                                                                                                      |
| `cfdna transitions`      | Extract nth-order transition probabilities in specifiable parts of the fragments                                                                                                                                       |
| `cfdna wps-peaks`        | Estimate nucleosome peaks from windowed protection scores                                                                                                                                                              |
| `cfdna wps`              | Calculate windowed protection scores per position (independent of `wps-peaks`)                                                                                                                                         |


# Optional preparation of intervals to count midpoints in
cfdna prepare_windows ...



# Transition probabilities in the first 10bp from each end
cfdna transitions --bam $BAM --output-dir $OUT --orders 1 --frame nearest --positions '..10' --min-fragment-length $MINLENGTH --max-fragment-length $MAXLENGTH --gc $OUT/gc_bias --scaling-factors $OUT/coverage_weights/<prefix>.scaling_factors.tsv --blacklist $BLACKLIST --n-threads $THREADS

# End motifs
cfdna ends ... # NOT IMPLEMENTED YET (breakpoint motif example also?)

# Nucleosome peaks from windowed protection scores
cfdna wps-peaks --bam $BAM --output-dir $OUT/wps_peaks --min-fragment-length 120 --max-fragment-length 180 --window-size 120 --gc $OUT/gc_bias --scaling-factors $OUT/coverage_weights/<prefix>.scaling_factors.tsv --blacklist $BLACKLIST --n-threads $THREADS

# Statistics on nucleosome peaks per 5Mb
cfdna wps-peaks --bam $BAM --output-dir $OUT/wps_peaks_statistics_per_5mb --by-size 5000000 --per-window stats --min-fragment-length 120 --max-fragment-length 180 --window-size 120 --gc $OUT/gc_bias --scaling-factors $OUT/coverage_weights/<prefix>.scaling_factors.tsv --blacklist $BLACKLIST --n-threads $THREADS