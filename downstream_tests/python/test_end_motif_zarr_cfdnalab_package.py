from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cfdnalab


def test_cfdnalab_package_reads_dense_global_end_motifs(
    dense_global_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(dense_global_end_zarr_path)

    assert isinstance(end_motifs, cfdnalab.GlobalEndMotifCounts)
    assert end_motifs.storage_mode() == "dense"
    assert end_motifs.row_mode() == "global"
    assert end_motifs.motifs() == ["_A", "_C", "_G", "_T"]
    assert end_motifs.dense_counts_zarr_array().shape == (1, 4)
    np.testing.assert_allclose(
        end_motifs.dense_counts_vec(),
        np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        end_motifs.dense_data_frame(),
        pd.DataFrame(
            {
                "motif_index": np.array([0, 1, 2, 3], dtype=np.int32),
                "motif": np.array(["_A", "_C", "_G", "_T"], dtype=str),
                "count": np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float64),
            }
        ),
    )
    with pytest.raises(ValueError, match="only available for sparse_coo output"):
        end_motifs.sparse_coo_data_frame()


def test_cfdnalab_package_reads_sparse_windowed_end_motifs(
    sparse_windowed_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(sparse_windowed_end_zarr_path)

    assert isinstance(end_motifs, cfdnalab.WindowedEndMotifCounts)
    assert end_motifs.storage_mode() == "sparse_coo"
    assert end_motifs.row_mode() == "bed"
    assert end_motifs.motifs() == ["_A", "_G"]
    np.testing.assert_allclose(
        end_motifs.dense_counts_matrix(),
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_for_window(1),
        np.array([1.0, 0.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_for_motif("_G"),
        np.array([1.0, 0.0], dtype=np.float64),
    )

    window_coo = end_motifs.sparse_coo_for_window(1)
    assert window_coo.shape == (1, 2)
    np.testing.assert_array_equal(window_coo.row, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(window_coo.col, np.array([0], dtype=np.int32))
    np.testing.assert_allclose(window_coo.data, np.array([1.0], dtype=np.float64))

    windows = end_motifs.windows()
    assert windows["window_idx"].tolist() == [0, 1]
    assert windows["chromosome_name"].tolist() == ["chr1", "chr1"]
    assert windows["window_start_bp"].tolist() == [10, 19]
    assert windows["window_end_bp"].tolist() == [11, 20]
    assert windows["blacklisted_fraction"].tolist() == [0.0, 0.0]


def test_cfdnalab_package_reads_sparse_grouped_end_motifs(
    sparse_grouped_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(sparse_grouped_end_zarr_path)

    assert isinstance(end_motifs, cfdnalab.GroupedEndMotifCounts)
    assert end_motifs.storage_mode() == "sparse_coo"
    assert end_motifs.row_mode() == "grouped_bed"
    assert end_motifs.motifs() == ["_A", "_G"]
    assert end_motifs.group_idx("alpha") == 1
    np.testing.assert_allclose(
        end_motifs.dense_counts_matrix(),
        np.array([[1.0, 2.0], [1.0, 0.0], [0.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_for_group("beta"),
        np.array([1.0, 2.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_for_motif("_A"),
        np.array([1.0, 1.0, 0.0], dtype=np.float64),
    )

    pd.testing.assert_frame_equal(
        end_motifs.groups(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1, 2], dtype=np.int32),
                "group_name": np.array(["beta", "alpha", "gamma"], dtype=str),
                "eligible_windows": np.array([2, 1, 1], dtype=np.uint32),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.0], dtype=np.float32),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        end_motifs.sparse_coo_data_frame(),
        pd.DataFrame(
            {
                "row": np.array([0, 0, 1], dtype=np.uint64),
                "motif_index": np.array([0, 1, 0], dtype=np.uint64),
                "motif": np.array(["_A", "_G", "_A"], dtype=str),
                "count": np.array([1.0, 2.0, 1.0], dtype=np.float64),
            }
        ),
    )

    beta_coo = end_motifs.sparse_coo_for_group("beta")
    assert beta_coo.shape == (1, 2)
    np.testing.assert_array_equal(beta_coo.row, np.array([0, 0], dtype=np.int32))
    np.testing.assert_array_equal(beta_coo.col, np.array([0, 1], dtype=np.int32))
    np.testing.assert_allclose(beta_coo.data, np.array([1.0, 2.0], dtype=np.float64))

    alpha_frame = end_motifs.dense_data_frame_for_group("alpha")
    assert alpha_frame["group_name"].unique().tolist() == ["alpha"]
    assert alpha_frame["motif"].tolist() == ["_A", "_G"]
    assert alpha_frame["count"].tolist() == [1.0, 0.0]
