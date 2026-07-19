from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cfdnalab


def test_cfdnalab_package_reads_dense_global_ref_kmers(
    dense_global_ref_kmer_zarr_path: Path,
) -> None:
    ref_kmers = cfdnalab.read_ref_kmers(dense_global_ref_kmer_zarr_path)

    assert isinstance(ref_kmers, cfdnalab.GlobalRefKmerFrequencies)
    assert ref_kmers.storage_mode() == "dense"
    assert ref_kmers.row_mode() == "global"
    assert ref_kmers.motif_axis_kind() == "motif"
    assert ref_kmers.kmer_size() == 3
    assert ref_kmers.canonical()
    assert ref_kmers.orientation() == "both"
    assert ref_kmers.all_motifs()
    assert ref_kmers.assign_by() == "count-overlap"
    assert [entry["name"] for entry in ref_kmers.reference_contig_footprint()] == [
        "chr1",
        "chr2",
    ]
    assert len(ref_kmers.motifs_metadata()) == 32
    assert ref_kmers.motifs_metadata()["motif"].head(4).tolist() == [
        "AAA",
        "AAC",
        "AAG",
        "AAT",
    ]
    assert ref_kmers.motifs_metadata()["motif"].tail(4).tolist() == [
        "TCA",
        "TCC",
        "TCG",
        "TCT",
    ]
    np.testing.assert_allclose(
        ref_kmers.dense_frequencies_array(motifs=["AAA", "ACG", "TAC", "CAT"]),
        np.array([[4 / 36, 7 / 36, 6 / 36, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(motifs=["AAA", "ACG", "TAC", "CAT"]),
        np.array([[4.0, 7.0, 6.0, 0.0]], dtype=np.float64),
    )
    assert ref_kmers.row_scaling_factors()["row_scaling_factor"].tolist() == [36.0]
    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(motifs=["TAC", "AAA"]),
        pd.DataFrame(
            {
                "row_label": np.array(["global", "global"], dtype=object),
                "motif_index": np.array([25, 0], dtype=np.int32),
                "motif": np.array(["TAC", "AAA"], dtype=object),
                "frequency": np.array([6 / 36, 4 / 36], dtype=np.float64),
                "count": np.array([6.0, 4.0], dtype=np.float64),
            }
        ),
        check_exact=False,
        rtol=1e-8,
    )


def test_cfdnalab_package_reads_sparse_windowed_ref_kmers(
    sparse_windowed_ref_kmer_zarr_path: Path,
) -> None:
    ref_kmers = cfdnalab.read_ref_kmers(sparse_windowed_ref_kmer_zarr_path)

    assert isinstance(ref_kmers, cfdnalab.WindowedRefKmerFrequencies)
    assert ref_kmers.storage_mode() == "sparse_coo"
    assert ref_kmers.row_mode() == "bed"
    assert ref_kmers.kmer_size() == 3
    assert ref_kmers.orientation() == "both"
    assert ref_kmers.motifs_metadata()["motif"].tolist() == [
        "CGT",
        "AAA",
        "TAC",
        "CCC",
        "GGG",
        "TTT",
        "ACG",
        "GTA",
    ]
    assert ref_kmers.has_motif("TTT")
    with pytest.raises(ValueError, match="would turn sparse reference k-mer output"):
        ref_kmers.dense_frequencies_array()
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(allow_densify=True),
        np.array(
            [
                [0.0, 1 / 2, 0.0, 1.0, 1.0, 1 / 2, 0.0, 0.0],
                [1.0, 0.0, 1.0, 1 / 6, 1 / 6, 0.0, 1.0, 1.0],
                [1 / 3, 0.0, 1 / 6, 2 / 3, 2 / 3, 0.0, 1 / 3, 1 / 6],
                [4 / 3, 0.0, 3 / 2, 0.0, 0.0, 0.0, 4 / 3, 3 / 2],
            ],
            dtype=np.float64,
        ),
    )
    np.testing.assert_allclose(
        ref_kmers.row_scaling_factors()["row_scaling_factor"].to_numpy(),
        np.array([3.0, 13 / 3, 7 / 3, 17 / 3], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.window_metadata(),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1, 2, 3], dtype=np.int32),
                "chrom": np.array(["chr1", "chr1", "chr2", "chr2"], dtype=object),
                "start": np.array([0, 8, 2, 12], dtype=np.int64),
                "end": np.array([9, 16, 13, 20], dtype=np.int64),
                "blacklisted_fraction": np.array(
                    [0.0, 1 / 8, 2 / 11, 0.0], dtype=np.float64
                ),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(
            window_idxs=[3, 0],
            motifs=["GTA", "AAA"],
            densify=True,
        ),
        pd.DataFrame(
            {
                "window_idx": np.array([3, 3, 0, 0], dtype=np.int32),
                "chrom": np.array(["chr2", "chr2", "chr1", "chr1"], dtype=object),
                "start": np.array([12, 12, 0, 0], dtype=np.int64),
                "end": np.array([20, 20, 9, 9], dtype=np.int64),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.0, 0.0]),
                "motif_index": np.array([7, 1, 7, 1], dtype=np.int32),
                "motif": np.array(["GTA", "AAA", "GTA", "AAA"], dtype=object),
                "frequency": np.array([9 / 34, 0.0, 0.0, 1 / 6], dtype=np.float64),
                "count": np.array([3 / 2, 0.0, 0.0, 1 / 2], dtype=np.float64),
            }
        ),
        check_exact=False,
        rtol=1e-8,
    )


def test_cfdnalab_package_reads_sparse_grouped_ref_kmers(
    sparse_grouped_ref_kmer_zarr_path: Path,
) -> None:
    ref_kmers = cfdnalab.read_ref_kmers(sparse_grouped_ref_kmer_zarr_path)

    assert isinstance(ref_kmers, cfdnalab.GroupedRefKmerFrequencies)
    assert ref_kmers.storage_mode() == "sparse_coo"
    assert ref_kmers.row_mode() == "grouped_bed"
    assert ref_kmers.orientation() == "both"
    assert ref_kmers.group_idx("alpha") == 1
    assert ref_kmers.motifs_metadata()["motif"].tolist() == [
        "CGT",
        "AAA",
        "TAC",
        "CCC",
        "GGG",
        "TTT",
        "ACG",
        "GTA",
    ]
    pd.testing.assert_frame_equal(
        ref_kmers.group_metadata(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1], dtype=np.int32),
                "group_name": np.array(["beta", "alpha"], dtype=object),
                "eligible_windows": np.array([2, 2], dtype=np.int32),
                "blacklisted_fraction": np.array([0.1, 1 / 16], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_allclose(
        ref_kmers.sparse_counts_matrix(
            groups=["alpha", "beta"],
            motifs=["GTA", "AAA"],
        ).toarray(),
        np.array([[5 / 2, 0.0], [1 / 6, 1 / 2]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(),
        pd.DataFrame(
            {
                "group_idx": np.array([0] * 8 + [1] * 6, dtype=np.int32),
                "group_name": np.array(["beta"] * 8 + ["alpha"] * 6, dtype=object),
                "eligible_windows": np.array([2] * 14, dtype=np.int32),
                "blacklisted_fraction": np.array([0.1] * 8 + [1 / 16] * 6),
                "motif_index": np.array(
                    [0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 3, 4, 6, 7],
                    dtype=np.int32,
                ),
                "motif": np.array(
                    [
                        "CGT",
                        "AAA",
                        "TAC",
                        "CCC",
                        "GGG",
                        "TTT",
                        "ACG",
                        "GTA",
                        "CGT",
                        "TAC",
                        "CCC",
                        "GGG",
                        "ACG",
                        "GTA",
                    ],
                    dtype=object,
                ),
                "frequency": np.array(
                    [
                        1 / 16,
                        3 / 32,
                        1 / 32,
                        5 / 16,
                        5 / 16,
                        3 / 32,
                        1 / 16,
                        1 / 32,
                        7 / 30,
                        1 / 4,
                        1 / 60,
                        1 / 60,
                        7 / 30,
                        1 / 4,
                    ],
                    dtype=np.float64,
                ),
                "count": np.array(
                    [
                        1 / 3,
                        1 / 2,
                        1 / 6,
                        5 / 3,
                        5 / 3,
                        1 / 2,
                        1 / 3,
                        1 / 6,
                        7 / 3,
                        5 / 2,
                        1 / 6,
                        1 / 6,
                        7 / 3,
                        5 / 2,
                    ],
                    dtype=np.float64,
                ),
            }
        ),
        check_exact=False,
        rtol=1e-8,
    )


def test_cfdnalab_package_reads_dense_grouped_motif_group_ref_kmers(
    dense_grouped_motif_group_ref_kmer_zarr_path: Path,
) -> None:
    ref_kmers = cfdnalab.read_ref_kmers(dense_grouped_motif_group_ref_kmer_zarr_path)

    assert isinstance(ref_kmers, cfdnalab.GroupedRefKmerFrequencies)
    assert ref_kmers.storage_mode() == "dense"
    assert ref_kmers.row_mode() == "grouped_bed"
    assert ref_kmers.motif_axis_kind() == "motif_group"
    assert ref_kmers.orientation() == "both"
    assert ref_kmers.all_motifs()
    assert ref_kmers.motifs_metadata()["motif"].tolist() == [
        "absent",
        "edge",
        "gc_rich",
        "homopolymer",
        "transition",
    ]
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(),
        np.array(
            [
                [0.0, 1 / 3, 10 / 3, 1.0, 2 / 3],
                [0.0, 5.0, 1 / 3, 0.0, 14 / 3],
            ],
            dtype=np.float64,
        ),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(groups=["alpha", "beta"], motifs=["edge", "absent"]),
        np.array([[5.0, 0.0], [1 / 3, 0.0]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(groups=["alpha", "beta"], motifs=["edge", "absent"]),
        pd.DataFrame(
            {
                "group_idx": np.array([1, 1, 0, 0], dtype=np.int32),
                "group_name": np.array(["alpha", "alpha", "beta", "beta"], dtype=object),
                "eligible_windows": np.array([2, 2, 2, 2], dtype=np.int32),
                "blacklisted_fraction": np.array([1 / 16, 1 / 16, 0.1, 0.1], dtype=np.float64),
                "motif_index": np.array([1, 0, 1, 0], dtype=np.int32),
                "motif": np.array(["edge", "absent", "edge", "absent"], dtype=object),
                "frequency": np.array([1 / 2, 0.0, 1 / 16, 0.0], dtype=np.float64),
                "count": np.array([5.0, 0.0, 1 / 3, 0.0], dtype=np.float64),
            }
        ),
        check_exact=False,
        rtol=1e-8,
    )
