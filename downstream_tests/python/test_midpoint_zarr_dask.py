from __future__ import annotations

from pathlib import Path

import dask.array as da
import numpy as np


def test_dask_reads_count_tensor_chunks(midpoint_zarr_path: Path) -> None:
    counts = da.from_zarr(str(midpoint_zarr_path), component="counts")

    assert counts.shape == (3, 3, 5)
    assert counts.dtype == np.dtype("float32")
    np.testing.assert_allclose(
        counts[1, :, :].compute(),
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [1.5, 0.5, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.5, 0.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )
