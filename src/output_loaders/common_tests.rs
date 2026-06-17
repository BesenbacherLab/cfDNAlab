use super::{DenseMatrix, contiguous_index_span};

/// Verify that dense matrix row iteration follows row-major matrix order.
#[test]
fn dense_matrix_rows_iterates_rows_in_matrix_order() -> anyhow::Result<()> {
    // Arrange: row-major values represent two rows with three columns each.
    let matrix = DenseMatrix::from_row_major(vec![1, 2, 3, 4, 5, 6], 2, 3)?;

    // Act.
    let rows = matrix
        .rows()
        .map(|row| row.to_vec())
        .collect::<Vec<_>>();

    // Assert.
    assert_eq!(rows, vec![vec![1, 2, 3], vec![4, 5, 6]]);
    Ok(())
}

/// Verify that dense matrix row iteration preserves zero-column rows.
#[test]
fn dense_matrix_rows_handles_zero_column_rows() -> anyhow::Result<()> {
    // Arrange:
    // A matrix can have rows with no columns. Iterating rows must still produce
    // one empty slice per row and stop at the declared row count.
    let matrix = DenseMatrix::<u32>::from_row_major(Vec::new(), 2, 0)?;

    // Act.
    let rows = matrix.rows().collect::<Vec<_>>();

    // Assert.
    assert_eq!(matrix.row(0), Some(&[][..]));
    assert_eq!(matrix.row(1), Some(&[][..]));
    assert_eq!(matrix.row(2), None);
    assert_eq!(rows, vec![&[] as &[u32], &[] as &[u32]]);
    Ok(())
}

/// Verify that an empty index selection has no contiguous span.
#[test]
fn contiguous_index_span_returns_none_for_empty_input() {
    // Arrange: no first index exists, so there is no half-open range to return.
    let indices = [];

    // Act.
    let span = contiguous_index_span(&indices);

    // Assert.
    assert_eq!(span, None);
}

/// Verify that a single selected index maps to a one-column span.
#[test]
fn contiguous_index_span_returns_single_index_range() {
    // Arrange: one selected index is contiguous with itself.
    let indices = [4];

    // Act.
    let span = contiguous_index_span(&indices);

    // Assert: index 4 alone covers the half-open column range 4..5.
    assert_eq!(span, Some((4, 5)));
}

/// Verify that ascending adjacent indices map to one half-open span.
#[test]
fn contiguous_index_span_returns_half_open_range_for_contiguous_indices() {
    // Arrange: selected indices 2, 3, and 4 form one contiguous run.
    let indices = [2, 3, 4];

    // Act.
    let span = contiguous_index_span(&indices);

    // Assert: the run starts at 2 and ends just after the final selected index.
    assert_eq!(span, Some((2, 5)));
}

/// Verify that gaps prevent a selection from being represented as one span.
#[test]
fn contiguous_index_span_rejects_gapped_indices() {
    // Arrange: index 3 is missing, so this cannot be copied as one range.
    let indices = [2, 4, 5];

    // Act.
    let span = contiguous_index_span(&indices);

    // Assert.
    assert_eq!(span, None);
}

/// Verify that duplicate indices are not treated as contiguous.
#[test]
fn contiguous_index_span_rejects_duplicate_indices() {
    // Arrange: the second selected value repeats index 2 instead of advancing to 3.
    let indices = [2, 2, 3];

    // Act.
    let span = contiguous_index_span(&indices);

    // Assert.
    assert_eq!(span, None);
}

/// Verify that decreasing indices are not treated as one contiguous span.
#[test]
fn contiguous_index_span_rejects_decreasing_indices() {
    // Arrange: the selected indices are adjacent as a set, but not in ascending order.
    let indices = [4, 3, 2];

    // Act.
    let span = contiguous_index_span(&indices);

    // Assert.
    assert_eq!(span, None);
}
