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
    pd.testing.assert_frame_equal(
        end_motifs.motifs_metadata(),
        pd.DataFrame(
            {
                "motif_index": np.array([0, 1, 2, 3], dtype=np.int32),
                "motif": np.array(["_A", "_C", "_G", "_T"], dtype=object),
            }
        ),
    )
    assert end_motifs.dense_counts_zarr_array().shape == (1, 4)
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(),
        np.array([[1.0, 0.0, 1.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(motifs=["_G", "_A"]),
        np.array([[1.0, 1.0]], dtype=np.float64),
    )
    assert end_motifs.data_frame(motifs=["_G", "_A"])["motif"].tolist() == ["_G", "_A"]
    pd.testing.assert_frame_equal(
        end_motifs.data_frame(),
        pd.DataFrame(
            {
                "row_label": np.array(
                    ["global", "global", "global", "global"], dtype=object
                ),
                "motif_index": np.array([0, 1, 2, 3], dtype=np.int32),
                "motif": np.array(["_A", "_C", "_G", "_T"], dtype=str),
                "count": np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float64),
            }
        ),
    )


def test_cfdnalab_package_reads_sparse_windowed_end_motifs(
    sparse_windowed_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(sparse_windowed_end_zarr_path)

    assert isinstance(end_motifs, cfdnalab.WindowedEndMotifCounts)
    assert end_motifs.storage_mode() == "sparse_coo"
    assert end_motifs.row_mode() == "bed"
    assert end_motifs.motifs_metadata()["motif"].tolist() == ["_A", "_G"]
    with pytest.raises(ValueError, match="would densify a sparse end-motif store"):
        end_motifs.dense_counts_array()
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(allow_densify=True),
        np.array([[0.0, 1.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(window_idxs=1, allow_densify=True),
        np.array([[1.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(motifs="_G", allow_densify=True),
        np.array([[1.0], [0.0], [1.0]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        end_motifs.data_frame(
            motifs="_A",
            densify=True,
            max_blacklisted_fraction=0.0,
        ),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1, 2], dtype=np.int32),
                "chrom": np.array(["chr1", "chr1", "chr2"], dtype=object),
                "start": np.array([10, 19, 10], dtype=np.int64),
                "end": np.array([11, 20, 11], dtype=np.int64),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.0], dtype=np.float64),
                "motif_index": np.array([0, 0, 0], dtype=np.int32),
                "motif": ["_A", "_A", "_A"],
                "count": np.array([0.0, 1.0, 0.0], dtype=np.float64),
            }
        ),
    )

    window_coo = end_motifs.sparse_counts_matrix(window_idxs=1)
    assert window_coo.shape == (1, 2)
    np.testing.assert_array_equal(window_coo.row, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(window_coo.col, np.array([0], dtype=np.int32))
    np.testing.assert_allclose(window_coo.data, np.array([1.0], dtype=np.float64))

    motif_coo = end_motifs.sparse_counts_matrix(motifs="_G")
    assert motif_coo.shape == (3, 1)
    np.testing.assert_array_equal(motif_coo.row, np.array([0, 2], dtype=np.int32))
    np.testing.assert_array_equal(motif_coo.col, np.array([0, 0], dtype=np.int32))
    np.testing.assert_allclose(motif_coo.data, np.array([1.0, 1.0], dtype=np.float64))
    np.testing.assert_allclose(
        end_motifs.sparse_counts_matrix(
            window_idxs=[1, 0],
            motifs=["_G", "_A"],
        ).toarray(),
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.float64),
    )
    ordered_dense = end_motifs.data_frame(
        window_idxs=[1, 0],
        motifs=["_G", "_A"],
        densify=True,
    )
    assert ordered_dense["window_idx"].tolist() == [1, 1, 0, 0]
    assert ordered_dense["motif"].tolist() == ["_G", "_A", "_G", "_A"]
    assert ordered_dense["count"].tolist() == [0.0, 1.0, 1.0, 0.0]

    windows = end_motifs.window_metadata()
    assert windows["window_idx"].tolist() == [0, 1, 2]
    assert windows["chrom"].tolist() == ["chr1", "chr1", "chr2"]
    assert windows["start"].tolist() == [10, 19, 10]
    assert windows["end"].tolist() == [11, 20, 11]
    assert windows["blacklisted_fraction"].tolist() == [0.0, 0.0, 0.0]


def test_cfdnalab_package_reads_sparse_windowed_selected_motif_file_end_motifs(
    sparse_windowed_selected_motifs_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(
        sparse_windowed_selected_motifs_end_zarr_path
    )

    assert isinstance(end_motifs, cfdnalab.WindowedEndMotifCounts)
    assert end_motifs.storage_mode() == "sparse_coo"
    assert end_motifs.row_mode() == "bed"
    assert end_motifs.motifs_metadata()["motif"].tolist() == [
        "GT_AC",
        "AC_GT",
    ]
    assert not end_motifs.has_motif("TT_TT")
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(allow_densify=True),
        np.array(
            [[0.0, 1.0], [1.0, 0.0], [0.0, 1.0]],
            dtype=np.float64,
        ),
    )
    np.testing.assert_allclose(
        end_motifs.sparse_counts_matrix(motifs=["AC_GT", "GT_AC"]).toarray(),
        np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]], dtype=np.float64),
    )
    with pytest.raises(KeyError, match="Unknown end-motif label"):
        end_motifs.data_frame(motifs="TT_TT", densify=True)


def test_cfdnalab_package_reads_sparse_grouped_end_motifs(
    sparse_grouped_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(sparse_grouped_end_zarr_path)

    assert isinstance(end_motifs, cfdnalab.GroupedEndMotifCounts)
    assert end_motifs.storage_mode() == "sparse_coo"
    assert end_motifs.row_mode() == "grouped_bed"
    assert end_motifs.motifs_metadata()["motif"].tolist() == ["_A", "_G"]
    assert end_motifs.group_idx("alpha") == 1
    with pytest.raises(ValueError, match="would densify a sparse end-motif store"):
        end_motifs.dense_counts_array()
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(allow_densify=True),
        np.array([[1.0, 2.0], [1.0, 0.0], [0.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(groups="beta", allow_densify=True),
        np.array([[1.0, 2.0]], dtype=np.float64),
    )
    assert end_motifs.data_frame(
        groups="beta",
        densify=True,
        max_blacklisted_fraction=0.0,
    )["count"].tolist() == [1.0, 2.0]
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(motifs="_A", allow_densify=True),
        np.array([[1.0], [1.0], [0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.sparse_counts_matrix(
            groups=["alpha", "beta"],
            motifs=["_G", "_A"],
        ).toarray(),
        np.array([[0.0, 1.0], [2.0, 1.0]], dtype=np.float64),
    )

    pd.testing.assert_frame_equal(
        end_motifs.group_metadata(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1, 2], dtype=np.int32),
                "group_name": np.array(["beta", "alpha", "gamma"], dtype=str),
                "eligible_windows": np.array([2, 1, 1], dtype=np.int32),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.0], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        end_motifs.data_frame(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 0, 1], dtype=np.int32),
                "group_name": np.array(["beta", "beta", "alpha"], dtype=str),
                "eligible_windows": np.array([2, 2, 1], dtype=np.int32),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.0], dtype=np.float64),
                "motif_index": np.array([0, 1, 0], dtype=np.int32),
                "motif": np.array(["_A", "_G", "_A"], dtype=str),
                "count": np.array([1.0, 2.0, 1.0], dtype=np.float64),
            }
        ),
    )

    beta_coo = end_motifs.sparse_counts_matrix(groups="beta")
    assert beta_coo.shape == (1, 2)
    np.testing.assert_array_equal(beta_coo.row, np.array([0, 0], dtype=np.int32))
    np.testing.assert_array_equal(beta_coo.col, np.array([0, 1], dtype=np.int32))
    np.testing.assert_allclose(beta_coo.data, np.array([1.0, 2.0], dtype=np.float64))

    alpha_frame = end_motifs.data_frame(groups="alpha", densify=True)
    assert alpha_frame["group_name"].unique().tolist() == ["alpha"]
    assert alpha_frame["motif"].tolist() == ["_A", "_G"]
    assert alpha_frame["count"].tolist() == [1.0, 0.0]

    assert end_motifs.data_frame(groups="gamma").empty
    gamma_dense = end_motifs.data_frame(groups="gamma", densify=True)
    assert gamma_dense["group_name"].tolist() == ["gamma", "gamma"]
    assert gamma_dense["motif"].tolist() == ["_A", "_G"]
    assert gamma_dense["count"].tolist() == [0.0, 0.0]


def test_cfdnalab_package_reads_sparse_grouped_motif_group_end_motifs(
    sparse_grouped_motif_group_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(sparse_grouped_motif_group_end_zarr_path)

    assert isinstance(end_motifs, cfdnalab.GroupedEndMotifCounts)
    assert end_motifs.storage_mode() == "sparse_coo"
    assert end_motifs.row_mode() == "grouped_bed"
    pd.testing.assert_frame_equal(
        end_motifs.motifs_metadata(),
        pd.DataFrame(
            {
                "motif_index": np.array([0, 1], dtype=np.int32),
                "motif": np.array(["left-hit", "right-hit"], dtype=object),
            }
        ),
    )
    assert end_motifs.motif_idx("right-hit") == 1
    assert end_motifs.has_motif("left-hit")
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(allow_densify=True),
        np.array([[2.0, 1.0], [0.0, 1.0], [0.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.sparse_counts_matrix(
            groups=["alpha", "beta"],
            motifs=["right-hit", "left-hit"],
        ).toarray(),
        np.array([[1.0, 0.0], [1.0, 2.0]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        end_motifs.data_frame(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 0, 1], dtype=np.int32),
                "group_name": np.array(["beta", "beta", "alpha"], dtype=str),
                "eligible_windows": np.array([2, 2, 1], dtype=np.int32),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.0], dtype=np.float64),
                "motif_index": np.array([0, 1, 1], dtype=np.int32),
                "motif": np.array(
                    ["left-hit", "right-hit", "right-hit"],
                    dtype=object,
                ),
                "count": np.array([2.0, 1.0, 1.0], dtype=np.float64),
            }
        ),
    )
    alpha_frame = end_motifs.data_frame(groups="alpha", densify=True)
    assert alpha_frame["motif"].tolist() == ["left-hit", "right-hit"]
    assert alpha_frame["count"].tolist() == [0.0, 1.0]
    gamma_frame = end_motifs.data_frame(groups="gamma", densify=True)
    assert gamma_frame["motif"].tolist() == ["left-hit", "right-hit"]
    assert gamma_frame["count"].tolist() == [0.0, 0.0]
    with pytest.raises(KeyError, match="Unknown end-motif label"):
        end_motifs.data_frame(motifs="_A")


def test_cfdnalab_package_reads_sparse_grouped_wide_motif_group_end_motifs(
    sparse_grouped_wide_motif_group_end_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(
        sparse_grouped_wide_motif_group_end_zarr_path
    )

    assert isinstance(end_motifs, cfdnalab.GroupedEndMotifCounts)
    assert end_motifs.storage_mode() == "sparse_coo"
    assert end_motifs.row_mode() == "grouped_bed"
    pd.testing.assert_frame_equal(
        end_motifs.motifs_metadata(),
        pd.DataFrame(
            {
                "motif_index": np.array([0, 1], dtype=np.int32),
                "motif": np.array(["right-hit-wide", "left-hit-wide"], dtype=object),
            }
        ),
    )
    assert end_motifs.motif_idx("left-hit-wide") == 1
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(allow_densify=True),
        np.array([[1.0, 2.0], [1.0, 0.0], [0.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.sparse_counts_matrix(
            groups=["alpha", "beta"],
            motifs=["left-hit-wide", "right-hit-wide"],
        ).toarray(),
        np.array([[0.0, 1.0], [2.0, 1.0]], dtype=np.float64),
    )
    beta_frame = end_motifs.data_frame(groups="beta", densify=True)
    assert beta_frame["motif"].tolist() == ["right-hit-wide", "left-hit-wide"]
    assert beta_frame["count"].tolist() == [1.0, 2.0]
    with pytest.raises(KeyError, match="Unknown end-motif label"):
        end_motifs.data_frame(motifs="GT_AC")
