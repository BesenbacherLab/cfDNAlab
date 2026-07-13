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
    with pytest.raises(ValueError, match="two_sided_correction requires ref_kmers"):
        end_motifs.data_frame(two_sided_correction="joint")


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


def test_cfdnalab_package_corrects_two_sided_end_motifs_without_same_motifs_file(
    sparse_windowed_two_sided_end_zarr_path: Path,
    sparse_windowed_end_motif_ref_kmer_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(sparse_windowed_two_sided_end_zarr_path)
    ref_kmers = cfdnalab.read_ref_kmers(sparse_windowed_end_motif_ref_kmer_zarr_path)

    assert end_motifs.motifs_metadata()["motif"].tolist() == ["AC_GT", "GT_AC"]
    assert ref_kmers.kmer_size() == 4
    assert ref_kmers.motifs_metadata()["motif"].tolist() == [
        "ACGT",
        "CGTA",
        "CGTT",
        "GTAC",
        "GTTT",
        "TACG",
        "TTTT",
    ]
    expected_counts = np.array(
        [[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]],
        dtype=np.float64,
    )
    # The stored columns are [AC_GT, GT_AC], and the three sample rows contain
    # [[1, 0], [0, 1], [1, 0]]. Seven positive reference k-mers make the joint
    # uniform frequency 1/7. Relative to uniform, ACGT frequency 1/4 gives
    # correction factor (1/4)/(1/7) = 7/4, while GTAC frequency 3/20 gives
    # (3/20)/(1/7) = 21/20. Dividing each observed count by its factor gives
    # [[4/7, 0], [0, 20/21], [4/7, 0]].
    expected_joint = np.array(
        [[4.0 / 7.0, 0.0], [0.0, 20.0 / 21.0], [4.0 / 7.0, 0.0]],
        dtype=np.float64,
    )
    # Five positive labels on each side make each side's uniform frequency
    # 1/5. Outside labels AC and GT both have frequency 1/4, so both factors
    # are (1/4)/(1/5) = 5/4. Inside labels GT and AC have frequencies 1/4 and
    # 3/20, giving factors 5/4 and 3/4. Split therefore divides AC_GT by
    # (5/4)*(5/4)=25/16 and GT_AC by (5/4)*(3/4)=15/16. Because each observed
    # count is 1, the corrected values are 16/25 and 16/15. Outside correction
    # divides by 5/4 for either label. Inside correction divides AC_GT by 5/4
    # and GT_AC by 3/4. In the stored row and column order, the split, outside,
    # and inside matrices are therefore [[16/25, 0], [0, 16/15], [16/25, 0]],
    # [[4/5, 0], [0, 4/5], [4/5, 0]], and
    # [[4/5, 0], [0, 4/3], [4/5, 0]], respectively.
    expected_split = np.array(
        [[16.0 / 25.0, 0.0], [0.0, 16.0 / 15.0], [16.0 / 25.0, 0.0]],
        dtype=np.float64,
    )
    expected_outside = np.array(
        [[4.0 / 5.0, 0.0], [0.0, 4.0 / 5.0], [4.0 / 5.0, 0.0]],
        dtype=np.float64,
    )
    expected_inside = np.array(
        [[4.0 / 5.0, 0.0], [0.0, 4.0 / 3.0], [4.0 / 5.0, 0.0]],
        dtype=np.float64,
    )

    with pytest.raises(ValueError, match="two-sided"):
        end_motifs.data_frame(ref_kmers=ref_kmers, densify=True)

    joint = end_motifs.data_frame(
        ref_kmers=ref_kmers,
        densify=True,
        two_sided_correction="joint",
    )
    assert joint["window_idx"].tolist() == [0, 0, 1, 1, 2, 2]
    assert joint["motif"].tolist() == ["AC_GT", "GT_AC"] * 3
    np.testing.assert_allclose(joint["count"].to_numpy(), expected_counts.ravel())
    np.testing.assert_allclose(
        joint["corrected_count"].to_numpy(),
        expected_joint.ravel(),
    )
    np.testing.assert_allclose(
        joint["corrected_frequency"].to_numpy(),
        expected_counts.ravel(),
    )

    np.testing.assert_allclose(
        end_motifs.corrected_counts_array(
            ref_kmers,
            allow_densify=True,
            two_sided_correction="split",
        ),
        expected_split,
    )
    np.testing.assert_allclose(
        end_motifs.sparse_corrected_counts_matrix(
            ref_kmers,
            two_sided_correction="split",
        ).toarray(),
        expected_split,
    )

    outside = end_motifs.data_frame(
        ref_kmers=ref_kmers,
        densify=True,
        two_sided_correction="outside",
    )
    assert outside["motif"].tolist() == ["AC_", "GT_"] * 3
    np.testing.assert_allclose(outside["count"].to_numpy(), expected_counts.ravel())
    np.testing.assert_allclose(
        outside["corrected_count"].to_numpy(),
        expected_outside.ravel(),
    )
    np.testing.assert_allclose(
        end_motifs.corrected_counts_array(
            ref_kmers,
            allow_densify=True,
            two_sided_correction="outside",
        ),
        expected_outside,
    )

    inside = end_motifs.data_frame(
        ref_kmers=ref_kmers,
        densify=True,
        two_sided_correction="inside",
    )
    assert inside["motif"].tolist() == ["_GT", "_AC"] * 3
    np.testing.assert_allclose(inside["count"].to_numpy(), expected_counts.ravel())
    np.testing.assert_allclose(
        inside["corrected_count"].to_numpy(),
        expected_inside.ravel(),
    )
    np.testing.assert_allclose(
        end_motifs.sparse_corrected_counts_matrix(
            ref_kmers,
            two_sided_correction="inside",
        ).toarray(),
        expected_inside,
    )


def test_cfdnalab_package_corrects_two_sided_end_motifs_with_same_motifs_file(
    sparse_windowed_selected_motifs_end_zarr_path: Path,
    sparse_windowed_selected_end_motifs_ref_kmer_zarr_path: Path,
) -> None:
    end_motifs = cfdnalab.read_end_motifs(
        sparse_windowed_selected_motifs_end_zarr_path
    )
    ref_kmers = cfdnalab.read_ref_kmers(
        sparse_windowed_selected_end_motifs_ref_kmer_zarr_path
    )

    assert end_motifs.motifs_metadata()["motif"].tolist() == ["GT_AC", "AC_GT"]
    assert ref_kmers.motifs_metadata()["motif"].tolist() == [
        "GTAC",
        "ACGT",
        "GTTT",
        "TTTT",
    ]
    expected_counts = np.array(
        [[0.0, 1.0], [1.0, 0.0], [0.0, 1.0]],
        dtype=np.float64,
    )
    # The stored columns are [GT_AC, AC_GT], and the three sample rows contain
    # [[0, 1], [1, 0], [0, 1]]. Four positive reference k-mers make the joint
    # uniform frequency 1/4. Relative to uniform, GTAC frequency 3/11 gives
    # correction factor (3/11)/(1/4) = 12/11, while ACGT frequency 5/11 gives
    # (5/11)/(1/4) = 20/11. Dividing each observed count by its factor gives
    # [[0, 11/20], [11/12, 0], [0, 11/20]].
    expected_joint = np.array(
        [[0.0, 11.0 / 20.0], [11.0 / 12.0, 0.0], [0.0, 11.0 / 20.0]],
        dtype=np.float64,
    )
    # Three positive labels on each side make each side's uniform frequency
    # 1/3. Outside labels AC and GT both have frequency 5/11, so both factors
    # are (5/11)/(1/3) = 15/11. Inside labels GT and AC have frequencies 5/11
    # and 3/11, giving factors 15/11 and 9/11. Split therefore divides AC_GT
    # by (15/11)*(15/11)=225/121 and GT_AC by
    # (15/11)*(9/11)=135/121. Because each observed count is 1, the corrected
    # values are 121/225 and 121/135. Outside correction divides by 15/11 for
    # either label. Inside correction divides AC_GT by 15/11 and GT_AC by 9/11.
    # In stored order, the split, outside, and inside matrices are therefore
    # [[0, 121/225], [121/135, 0], [0, 121/225]],
    # [[0, 11/15], [11/15, 0], [0, 11/15]], and
    # [[0, 11/15], [11/9, 0], [0, 11/15]], respectively.
    expected_split = np.array(
        [[0.0, 121.0 / 225.0], [121.0 / 135.0, 0.0], [0.0, 121.0 / 225.0]],
        dtype=np.float64,
    )
    expected_outside = np.array(
        [[0.0, 11.0 / 15.0], [11.0 / 15.0, 0.0], [0.0, 11.0 / 15.0]],
        dtype=np.float64,
    )
    expected_inside = np.array(
        [[0.0, 11.0 / 15.0], [11.0 / 9.0, 0.0], [0.0, 11.0 / 15.0]],
        dtype=np.float64,
    )

    joint = end_motifs.data_frame(
        ref_kmers=ref_kmers,
        densify=True,
        two_sided_correction="joint",
    )
    assert joint["motif"].tolist() == ["GT_AC", "AC_GT"] * 3
    np.testing.assert_allclose(joint["count"].to_numpy(), expected_counts.ravel())
    np.testing.assert_allclose(
        joint["corrected_count"].to_numpy(),
        expected_joint.ravel(),
    )
    np.testing.assert_allclose(
        joint["corrected_frequency"].to_numpy(),
        expected_counts.ravel(),
    )

    np.testing.assert_allclose(
        end_motifs.corrected_counts_array(
            ref_kmers,
            allow_densify=True,
            two_sided_correction="split",
        ),
        expected_split,
    )
    np.testing.assert_allclose(
        end_motifs.sparse_corrected_counts_matrix(
            ref_kmers,
            two_sided_correction="split",
        ).toarray(),
        expected_split,
    )

    outside = end_motifs.data_frame(
        ref_kmers=ref_kmers,
        densify=True,
        two_sided_correction="outside",
    )
    assert outside["motif"].tolist() == ["GT_", "AC_"] * 3
    np.testing.assert_allclose(outside["count"].to_numpy(), expected_counts.ravel())
    np.testing.assert_allclose(
        outside["corrected_count"].to_numpy(),
        expected_outside.ravel(),
    )
    np.testing.assert_allclose(
        end_motifs.corrected_counts_array(
            ref_kmers,
            allow_densify=True,
            two_sided_correction="outside",
        ),
        expected_outside,
    )

    selected_outside = end_motifs.data_frame(
        ref_kmers=ref_kmers,
        motifs="AC_",
        densify=True,
        two_sided_correction="outside",
    )
    assert selected_outside["motif"].tolist() == ["AC_", "AC_", "AC_"]
    np.testing.assert_allclose(
        selected_outside["corrected_count"].to_numpy(),
        np.array([11.0 / 15.0, 0.0, 11.0 / 15.0], dtype=np.float64),
    )

    inside = end_motifs.data_frame(
        ref_kmers=ref_kmers,
        densify=True,
        two_sided_correction="inside",
    )
    assert inside["motif"].tolist() == ["_AC", "_GT"] * 3
    np.testing.assert_allclose(inside["count"].to_numpy(), expected_counts.ravel())
    np.testing.assert_allclose(
        inside["corrected_count"].to_numpy(),
        expected_inside.ravel(),
    )
    np.testing.assert_allclose(
        end_motifs.sparse_corrected_counts_matrix(
            ref_kmers,
            two_sided_correction="inside",
        ).toarray(),
        expected_inside,
    )


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
                "motif": np.array(["left-hit-wide", "right-hit-wide"], dtype=object),
            }
        ),
    )
    assert end_motifs.motif_idx("left-hit-wide") == 0
    np.testing.assert_allclose(
        end_motifs.dense_counts_array(allow_densify=True),
        np.array([[2.0, 1.0], [0.0, 1.0], [0.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        end_motifs.sparse_counts_matrix(
            groups=["alpha", "beta"],
            motifs=["left-hit-wide", "right-hit-wide"],
        ).toarray(),
        np.array([[0.0, 1.0], [2.0, 1.0]], dtype=np.float64),
    )
    beta_frame = end_motifs.data_frame(groups="beta", densify=True)
    assert beta_frame["motif"].tolist() == ["left-hit-wide", "right-hit-wide"]
    assert beta_frame["count"].tolist() == [2.0, 1.0]
    with pytest.raises(KeyError, match="Unknown end-motif label"):
        end_motifs.data_frame(motifs="GT_AC")
