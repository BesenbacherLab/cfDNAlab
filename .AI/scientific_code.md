# Scientific Code

Use this file for reducers, aggregation, masking, tiling, normalization, and similar scientific or counting code.

Any reducer, aggregation, masking, tiling, or normalization change requires an invariant summary before edits.

Representation cleanups must not change persisted row identity or grouping keys.

If tests fail after recent edits, assume regression first and prove otherwise before claiming an old bug.

For scientific/counting code, prefer smaller patches and targeted verification over broad cleanup.
