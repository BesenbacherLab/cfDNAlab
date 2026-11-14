**Title:** `2_data_transformation.py --both`: right-end rows in `M__`/`T__` tables are dropped or padded with NaN

### Summary

With `-b/--both`, `BaseCompAbCalc` reuses a single `BCDf` data frame for both ends (`workflow/scripts/1_FrEIA/2_data_transformation.py`, lines 122–124 at `2c32fb9`):
- The first column assignment, `BCDf.CountLen = hBCDf.strand` (line 132), gives the empty frame the left end's index.
- Every later assignment, including all of the right end's, is aligned to that index by pandas.
- The right-end table therefore always has exactly as many rows as the left-end table, so:
  - **If the right end has more observed `(read_len, motif)` combinations**, its extra rows are silently dropped. These are the longest lengths, because `groupby` sorts by `read_len`. Its `CountLen` then no longer sums to the number of fragments.
  - **If it has fewer**, rows with `NaN` in `Base`, `Length`, `CountLen`, `CountSamp` and `RelAbSamp` are appended with `WhichEnd == "rightmost"`.

This affects the `M__*.pq` and `T__*.pq` outputs (`MonoNcAnalysis` and `TriNcAnalysis`) whenever `--both` is used. The Snakemake rule in `workflow/rules/2_FrEIA/2_FrEIA.smk` does not pass `--both`, so the default workflow outputs are unaffected.

### Minimal example

Three fragments, run through `BaseCompAbCalc(df, "g", ["leftmost_tnc_seq", "rightmost_tnc_seq"], "Genome_Lvl")`:

| read_len | leftmost_tnc_seq | rightmost_tnc_seq |
| --- | --- | --- |
| 150 | AAA | GGG |
| 150 | AAA | TTT |
| 160 | CCC | GGG |

The left end has 2 `(read_len, motif)` groups and the right end has 3. The right-end rows written:

| Length | Base | CountLen | Expected |
| --- | --- | --- | --- |
| 150 | GGG | 1 | 1 |
| 150 | TTT | 1 | 1 |
| *(missing)* | | | 160, GGG, 1 |

Swapping the two motif columns gives a third right-end row with `NaN` in every count column instead.

### Suggested fix

Create the frame inside the loop, so that each end gets its own index:

```python
def BaseCompAbCalc(SampleDf, group, whichEnd, lvl):
    OutDf = pd.DataFrame()

    for e in whichEnd:
        BCDf = pd.DataFrame(columns=["Unit", "Base", "Length", "CountSamp",
                                     "CountLen", "RelAbSamp", "WhichEnd",
                                     "WhichGroup"])
        ...
```

Found at commit `2c32fb9`; `workflow/envs/FrEIA_env.yaml` pins pandas 1.1.5. Reproduced with pandas 1.5.3. Both versions align a Series to the frame's existing index on column assignment.
