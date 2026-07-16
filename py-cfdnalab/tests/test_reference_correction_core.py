from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

import cfdnalab.reference_correction as reference_correction

# Keep this fixture and its rational expectations identical to the Rust and R
# core correction tests. It is the language-parity specification for the math.


def _shared_end_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "row_label": ["global"] * 4,
            "motif_index": np.arange(4, dtype=np.int64),
            "motif": ["A_C", "A_G", "T_C", "T_G"],
            "count": [2.0, 4.0, 6.0, 8.0],
        }
    )


def _shared_reference_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "row_label": ["global"] * 4,
            "reference_motif": ["AC", "AG", "TC", "TG"],
            "reference_frequency": [1.0 / 8.0, 1.0 / 8.0, 1.0 / 4.0, 1.0 / 2.0],
        }
    )


def _two_row_end_rows() -> pd.DataFrame:
    beta_rows = _shared_end_rows().assign(row_label="beta")
    alpha_rows = _shared_end_rows().assign(
        row_label="alpha",
        count=[2.0, 4.0, 0.0, 0.0],
    )
    return pd.concat([beta_rows, alpha_rows], ignore_index=True)


def _two_row_reference_rows() -> pd.DataFrame:
    alpha_rows = _shared_reference_rows().assign(
        row_label="alpha",
        reference_frequency=[0.5, 0.5, 0.0, 0.0],
    )
    beta_rows = _shared_reference_rows().assign(row_label="beta")
    return pd.concat([alpha_rows, beta_rows], ignore_index=True)


def _mode(
    mode: str,
    side_labels: list[str] | None = None,
) -> reference_correction._ReferenceCorrectionMode:
    return reference_correction._ReferenceCorrectionMode(
        mode=mode,
        outside_width=1,
        inside_width=1,
        side_labels=() if side_labels is None else tuple(side_labels),
    )


def _with_frequencies(corrected: pd.DataFrame) -> pd.DataFrame:
    return reference_correction._add_corrected_frequency(
        corrected,
        ["row_label"],
    )


def test_joint_core_uses_full_motif_frequencies() -> None:
    end_rows = _shared_end_rows()
    reference_rows = _shared_reference_rows()
    ends = SimpleNamespace(end_motifs=SimpleNamespace(motif_axis_kind="motif"))

    corrected = reference_correction._correct_exact_label_data_frame(
        ends,
        end_rows,
        reference_rows,
        ["row_label"],
        "error",
    )
    corrected = _with_frequencies(corrected)

    # Four positive reference motifs make the uniform frequency 1/4. Relative
    # to uniform, frequencies [1/8, 1/8, 1/4, 1/2] give correction factors
    # [1/2, 1/2, 1, 2] for [AC, AG, TC, TG]. Dividing original counts
    # [2, 4, 6, 8] by those factors gives [4, 8, 6, 4]. Their total is 22, so
    # dividing each corrected count by 22 gives [2/11, 4/11, 3/11, 2/11].
    np.testing.assert_allclose(
        corrected["reference_denominator"],
        [1.0 / 2.0, 1.0 / 2.0, 1.0, 2.0],
    )
    np.testing.assert_allclose(corrected["corrected_count"], [4.0, 8.0, 6.0, 4.0])
    np.testing.assert_allclose(
        corrected["corrected_frequency"],
        [2.0 / 11.0, 4.0 / 11.0, 3.0 / 11.0, 2.0 / 11.0],
    )


def test_joint_core_uses_support_from_each_reference_row() -> None:
    # Arrange: Sample rows are beta then alpha, while reference rows are alpha
    # then beta. Alpha supports two motifs and beta supports all four
    ends = SimpleNamespace(end_motifs=SimpleNamespace(motif_axis_kind="motif"))

    # Act
    corrected = reference_correction._correct_exact_label_data_frame(
        ends,
        _two_row_end_rows(),
        _two_row_reference_rows(),
        ["row_label"],
        "error",
    )

    # Assert: Beta's frequencies give denominators [1/2, 1/2, 1, 2]. Alpha's
    # two supported frequencies are both 1/2, giving [1, 1, 0, 0]. Output keeps
    # sample row and motif order, and unsupported zero counts remain zero
    assert corrected["row_label"].tolist() == ["beta"] * 4 + ["alpha"] * 4
    assert corrected["motif"].tolist() == ["A_C", "A_G", "T_C", "T_G"] * 2
    np.testing.assert_allclose(
        corrected["reference_denominator"],
        [0.5, 0.5, 1.0, 2.0, 1.0, 1.0, 0.0, 0.0],
    )
    np.testing.assert_allclose(
        corrected["corrected_count"],
        [4.0, 8.0, 6.0, 4.0, 2.0, 4.0, 0.0, 0.0],
    )


def test_split_core_multiplies_outside_and_inside_denominators() -> None:
    corrected = reference_correction._correct_split_data_frame(
        _shared_end_rows(),
        _shared_reference_rows(),
        ["row_label"],
        _mode("split"),
        "error",
    )
    corrected = _with_frequencies(corrected)

    # Two positive labels on each side make each side's uniform frequency 1/2.
    # Outside frequencies A=1/4 and T=3/4 give factors 1/2 and 3/2. Inside
    # frequencies C=3/8 and G=5/8 give factors 3/4 and 5/4. Multiplying matching
    # side factors gives [3/8, 5/8, 9/8, 15/8] for [A_C, A_G, T_C, T_G].
    # Dividing original counts [2, 4, 6, 8] by those factors gives
    # [16/3, 32/5, 16/3, 64/15]. Their total is 64/3, so normalization gives
    # frequencies [1/4, 3/10, 1/4, 1/5].
    np.testing.assert_allclose(
        corrected["reference_denominator"],
        [3.0 / 8.0, 5.0 / 8.0, 9.0 / 8.0, 15.0 / 8.0],
    )
    np.testing.assert_allclose(
        corrected["corrected_count"],
        [16.0 / 3.0, 32.0 / 5.0, 16.0 / 3.0, 64.0 / 15.0],
    )
    np.testing.assert_allclose(
        corrected["corrected_frequency"],
        [1.0 / 4.0, 3.0 / 10.0, 1.0 / 4.0, 1.0 / 5.0],
    )


def test_split_core_uses_side_support_from_each_reference_row() -> None:
    # Arrange: Alpha has only outside A support, while beta supports A and T.
    # Both rows support inside C and G. Reference and sample row order differ

    # Act
    corrected = reference_correction._correct_split_data_frame(
        _two_row_end_rows(),
        _two_row_reference_rows(),
        ["row_label"],
        _mode("split"),
        "error",
    )

    # Assert: Beta keeps the shared split denominators [3/8, 5/8, 9/8, 15/8].
    # Alpha has outside denominators A=1 and T=0, and inside denominators C=1
    # and G=1, giving full denominators [1, 1, 0, 0]
    assert corrected["row_label"].tolist() == ["beta"] * 4 + ["alpha"] * 4
    assert corrected["motif"].tolist() == ["A_C", "A_G", "T_C", "T_G"] * 2
    np.testing.assert_allclose(
        corrected["reference_denominator"],
        [3.0 / 8.0, 5.0 / 8.0, 9.0 / 8.0, 15.0 / 8.0, 1.0, 1.0, 0.0, 0.0],
    )
    np.testing.assert_allclose(
        corrected["corrected_count"],
        [16.0 / 3.0, 32.0 / 5.0, 16.0 / 3.0, 64.0 / 15.0, 2.0, 4.0, 0.0, 0.0],
    )


def test_split_core_handles_an_empty_sparse_reference_row() -> None:
    # Arrange: A sparse reference row with no stored motifs provides no outside
    # or inside support. Zero sample counts remain defined as zero
    end_rows = _shared_end_rows().assign(count=0.0)
    reference_rows = _shared_reference_rows().iloc[0:0]

    # Act
    corrected = reference_correction._correct_split_data_frame(
        end_rows,
        reference_rows,
        ["row_label"],
        _mode("split"),
        "error",
    )

    # Assert
    assert corrected["motif"].tolist() == ["A_C", "A_G", "T_C", "T_G"]
    np.testing.assert_allclose(corrected["reference_denominator"], np.zeros(4))
    np.testing.assert_allclose(corrected["corrected_count"], np.zeros(4))


def test_outside_core_aggregates_counts_before_correction() -> None:
    end_rows = _shared_end_rows()
    output_columns = end_rows.columns.tolist()
    end_rows["_cfdnalab_row_order"] = 0
    corrected = reference_correction._correct_side_data_frame(
        end_rows,
        _shared_reference_rows(),
        ["row_label"],
        _mode("outside", ["A_", "T_"]),
        output_columns,
        "error",
    )
    corrected = _with_frequencies(corrected)

    # Counts aggregate to A_=2+4=6 and T_=6+8=14. Two positive outside labels
    # make the uniform frequency 1/2. Relative to uniform, frequencies A=1/4
    # and T=3/4 give factors 1/2 and 3/2. Dividing the aggregated counts by
    # those factors gives [12, 28/3]. Their total is 64/3, so normalization
    # gives frequencies [9/16, 7/16].
    assert corrected["motif"].tolist() == ["A_", "T_"]
    np.testing.assert_allclose(corrected["count"], [6.0, 14.0])
    np.testing.assert_allclose(
        corrected["reference_denominator"],
        [1.0 / 2.0, 3.0 / 2.0],
    )
    np.testing.assert_allclose(corrected["corrected_count"], [12.0, 28.0 / 3.0])
    np.testing.assert_allclose(
        corrected["corrected_frequency"],
        [9.0 / 16.0, 7.0 / 16.0],
    )


def test_inside_core_aggregates_counts_before_correction() -> None:
    end_rows = _shared_end_rows()
    output_columns = end_rows.columns.tolist()
    end_rows["_cfdnalab_row_order"] = 0
    corrected = reference_correction._correct_side_data_frame(
        end_rows,
        _shared_reference_rows(),
        ["row_label"],
        _mode("inside", ["_C", "_G"]),
        output_columns,
        "error",
    )
    corrected = _with_frequencies(corrected)

    # Counts aggregate to _C=2+6=8 and _G=4+8=12. Two positive inside labels
    # make the uniform frequency 1/2. Relative to uniform, frequencies C=3/8
    # and G=5/8 give factors 3/4 and 5/4. Dividing the aggregated counts by
    # those factors gives [32/3, 48/5]. Their total is 304/15, so normalization
    # gives frequencies [10/19, 9/19].
    assert corrected["motif"].tolist() == ["_C", "_G"]
    np.testing.assert_allclose(corrected["count"], [8.0, 12.0])
    np.testing.assert_allclose(
        corrected["reference_denominator"],
        [3.0 / 4.0, 5.0 / 4.0],
    )
    np.testing.assert_allclose(corrected["corrected_count"], [32.0 / 3.0, 48.0 / 5.0])
    np.testing.assert_allclose(
        corrected["corrected_frequency"],
        [10.0 / 19.0, 9.0 / 19.0],
    )


def test_corrected_frequency_avoids_overflow_for_large_finite_counts() -> None:
    # Arrange: Adding two maximum finite floats overflows, but their equal
    # relative sizes still define frequencies of 1/2 and 1/2.
    maximum_float = np.finfo(np.float64).max
    corrected = pd.DataFrame(
        {
            "row_label": ["global", "global"],
            "corrected_count": [maximum_float, maximum_float],
        }
    )

    # Act
    corrected = reference_correction._add_corrected_frequency(
        corrected,
        ["row_label"],
    )

    # Assert
    np.testing.assert_allclose(corrected["corrected_frequency"], [0.5, 0.5])
