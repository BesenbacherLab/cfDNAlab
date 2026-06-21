use crate::interval::{IndexedInterval, Interval};
use fxhash::FxHashSet;

/// Reject public Zarr labels that cannot round-trip as one stable text value.
#[cfg(any(feature = "cmd_ends", feature = "cmd_midpoints"))]
pub(crate) fn validate_zarr_public_label(value: &str, field_name: &str) -> anyhow::Result<()> {
    anyhow::ensure!(
        !value.chars().any(char::is_control),
        "Zarr label {field_name} contains a control character"
    );
    Ok(())
}

/// A dense two-dimensional array stored in row-major order.
///
/// This stores a `Vec<T>` plus row and column counts. It keeps the output shape
/// next to the values without requiring a data frame or array crate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseMatrix<T> {
    values: Vec<T>,
    rows: usize,
    columns: usize,
}

impl<T> DenseMatrix<T> {
    /// Build a dense matrix from row-major values and explicit shape metadata.
    #[cfg(any(feature = "cmd_ends", feature = "cmd_lengths"))]
    pub(crate) fn from_row_major(
        values: Vec<T>,
        rows: usize,
        columns: usize,
    ) -> anyhow::Result<Self> {
        anyhow::ensure!(
            values.len() == rows.saturating_mul(columns),
            "dense matrix shape ({rows}, {columns}) does not match {} row-major values",
            values.len()
        );
        Ok(Self {
            values,
            rows,
            columns,
        })
    }

    /// Return the matrix shape as `(rows, columns)`.
    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.columns)
    }

    /// Return the number of matrix rows.
    pub fn row_count(&self) -> usize {
        self.rows
    }

    /// Return the number of matrix columns.
    pub fn column_count(&self) -> usize {
        self.columns
    }

    /// Return all matrix values in row-major order.
    pub fn values_row_major(&self) -> &[T] {
        &self.values
    }

    /// Return one matrix row, if `row_index` is in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index in the matrix.
    pub fn row(&self, row_index: usize) -> Option<&[T]> {
        if row_index >= self.rows {
            return None;
        }
        let start = row_index.checked_mul(self.columns)?;
        let end = start.checked_add(self.columns)?;
        self.values.get(start..end)
    }

    /// Return an iterator that yields one row of values at a time.
    ///
    /// Each item is a borrowed slice containing the values for one matrix row.
    /// The slice length is `column_count()`, rows are returned from first to
    /// last, and empty rows are returned as empty slices when the matrix has
    /// zero columns.
    ///
    /// ```no_run
    /// # use cfdnalab::output_loaders::load_lengths_output;
    /// # let lengths = load_lengths_output("sample.length_counts.tsv.zst")?;
    /// let row_totals = lengths
    ///     .counts()
    ///     .rows()
    ///     .map(|row_values| row_values.iter().copied().sum::<f64>())
    ///     .collect::<Vec<_>>();
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn rows(&self) -> impl ExactSizeIterator<Item = &[T]> + '_ {
        DenseMatrixRows {
            matrix: self,
            next_row: 0,
        }
    }

    /// Return one matrix value, if both indices are in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index in the matrix.
    /// - `column_index`:
    ///   Zero-based column index in the matrix.
    pub fn get(&self, row_index: usize, column_index: usize) -> Option<&T> {
        if row_index >= self.rows || column_index >= self.columns {
            return None;
        }
        let row_start = row_index.checked_mul(self.columns)?;
        let value_index = row_start.checked_add(column_index)?;
        self.values.get(value_index)
    }

    /// Return one row's values for a contiguous half-open column range.
    ///
    /// `start_column..end_column` uses zero-based column indices. Returns
    /// `None` if the row is out of bounds, if the column range is invalid, or
    /// if the range extends past the matrix column count.
    #[cfg(any(feature = "cmd_lengths", feature = "cmd_ends"))]
    pub(crate) fn row_values_for_column_range(
        &self,
        row_index: usize,
        start_column: usize,
        end_column: usize,
    ) -> Option<&[T]> {
        if row_index >= self.rows {
            return None;
        }
        if start_column > end_column || end_column > self.columns {
            return None;
        }
        let row_start = row_index.checked_mul(self.columns)?;
        let start = row_start.checked_add(start_column)?;
        let end = row_start.checked_add(end_column)?;
        self.values.get(start..end)
    }

    /// Consume the matrix and return the row-major value vector.
    pub fn into_values_row_major(self) -> Vec<T> {
        self.values
    }
}

/// Iterator returned by `DenseMatrix::rows()`.
///
/// The iterator borrows the matrix and tracks the next zero-based row index to
/// return. This keeps row iteration shape-aware, including matrices with rows
/// but no columns, where `slice::chunks()` cannot be used.
struct DenseMatrixRows<'a, T> {
    // Matrix whose row slices are returned
    matrix: &'a DenseMatrix<T>,
    // Zero-based row index returned by the next call to `Iterator::next`
    next_row: usize,
}

impl<'a, T> Iterator for DenseMatrixRows<'a, T> {
    /// Borrowed row slice yielded by the iterator.
    type Item = &'a [T];

    /// Return the next matrix row slice.
    fn next(&mut self) -> Option<Self::Item> {
        if self.next_row >= self.matrix.rows {
            return None;
        }
        let row_index = self.next_row;
        self.next_row += 1;
        self.matrix.row(row_index)
    }

    /// Return the exact number of row slices that remain.
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining_rows = self.len();
        (remaining_rows, Some(remaining_rows))
    }
}

impl<T> ExactSizeIterator for DenseMatrixRows<'_, T> {
    /// Return the number of row slices not yet yielded.
    fn len(&self) -> usize {
        self.matrix.rows.saturating_sub(self.next_row)
    }
}

/// One indexed half-open fragment length bin in output-column order.
pub type LengthBin = IndexedInterval<u32, usize>;

/// A dense three-dimensional array stored in row-major order.
///
/// This is used for output arrays where axes are part of the public command
/// result. Values are ordered with the last axis varying fastest.
#[cfg(feature = "cmd_midpoints")]
#[derive(Debug, Clone, PartialEq)]
pub struct DenseArray3<T> {
    values: Vec<T>,
    first: usize,
    second: usize,
    third: usize,
}

#[cfg(feature = "cmd_midpoints")]
impl<T> DenseArray3<T> {
    /// Build a dense three-dimensional array from row-major values and shape metadata.
    pub(crate) fn from_row_major(
        values: Vec<T>,
        first: usize,
        second: usize,
        third: usize,
    ) -> anyhow::Result<Self> {
        let expected_len = first
            .checked_mul(second)
            .and_then(|value| value.checked_mul(third))
            .ok_or_else(|| {
                anyhow::anyhow!("dense array shape ({first}, {second}, {third}) overflows usize")
            })?;
        anyhow::ensure!(
            values.len() == expected_len,
            "dense array shape ({first}, {second}, {third}) does not match {} row-major values",
            values.len()
        );
        Ok(Self {
            values,
            first,
            second,
            third,
        })
    }

    /// Return the array shape as `(first_axis, second_axis, third_axis)`.
    pub fn shape(&self) -> (usize, usize, usize) {
        (self.first, self.second, self.third)
    }

    /// Return all array values in row-major order.
    pub fn values_row_major(&self) -> &[T] {
        &self.values
    }

    /// Return one value, if all indices are in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `first_index`:
    ///   Zero-based index on the first axis.
    /// - `second_index`:
    ///   Zero-based index on the second axis.
    /// - `third_index`:
    ///   Zero-based index on the third axis.
    pub fn get(&self, first_index: usize, second_index: usize, third_index: usize) -> Option<&T> {
        if second_index >= self.second || third_index >= self.third {
            return None;
        }
        let first_offset = first_index.checked_mul(self.second)?;
        let second_offset = first_offset.checked_add(second_index)?;
        let row_start = second_offset.checked_mul(self.third)?;
        let value_index = row_start.checked_add(third_index)?;
        self.values.get(value_index)
    }

    /// Return all values along the third axis for fixed first and second indices.
    ///
    /// Parameters
    /// ----------
    /// - `first_index`:
    ///   Zero-based index on the first axis.
    /// - `second_index`:
    ///   Zero-based index on the second axis.
    pub fn values_along_third_axis(&self, first_index: usize, second_index: usize) -> Option<&[T]> {
        if second_index >= self.second {
            return None;
        }
        let first_offset = first_index.checked_mul(self.second)?;
        let second_offset = first_offset.checked_add(second_index)?;
        let start = second_offset.checked_mul(self.third)?;
        let end = start.checked_add(self.third)?;
        self.values.get(start..end)
    }

    /// Consume the array and return the row-major value vector.
    pub fn into_values_row_major(self) -> Vec<T> {
        self.values
    }
}

/// Metadata for one genomic window output row.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowRow {
    /// Zero-based row index in output-file order.
    pub index: usize,
    /// Chromosome or contig label from the output file.
    pub chrom: String,
    /// Checked half-open genomic interval for this row.
    pub interval: Interval<u64>,
    /// Fraction of the row interval covered by blacklisted bases, when written.
    pub blacklisted_fraction: Option<f64>,
}

/// Resolve an optional row selector to concrete zero-based source row indices.
pub(crate) fn resolve_row_indices(
    row_indices: Option<&[usize]>,
    row_count: usize,
    row_label: &str,
) -> anyhow::Result<Vec<usize>> {
    if let Some(row_indices) = row_indices {
        return Ok(row_indices.to_vec());
    }
    anyhow::ensure!(
        row_count > 0,
        "cannot select all {row_label} rows from an empty output"
    );
    Ok((0..row_count).collect())
}

/// Reject duplicate zero-based indices before a selection is copied.
pub(crate) fn ensure_unique_indices(indices: &[usize], label: &str) -> anyhow::Result<()> {
    let mut seen_indices = FxHashSet::default();
    for &index in indices {
        anyhow::ensure!(
            seen_indices.insert(index),
            "{label} indices contain duplicate value {index}"
        );
    }
    Ok(())
}

#[cfg(any(
    feature = "cmd_lengths",
    feature = "cmd_ends",
    feature = "cmd_midpoints"
))]
/// Reject duplicate labels before resolving label-based selectors.
pub(crate) fn ensure_unique_labels<S: AsRef<str>>(
    labels: &[S],
    label_name: &str,
) -> anyhow::Result<()> {
    let mut seen_labels = FxHashSet::default();
    for label in labels {
        let label = label.as_ref();
        anyhow::ensure!(
            seen_labels.insert(label),
            "{label_name} contains duplicate value '{label}'"
        );
    }
    Ok(())
}

#[cfg(any(
    feature = "cmd_lengths",
    feature = "cmd_ends",
    feature = "cmd_midpoints"
))]
/// Return the half-open span covered by strictly contiguous indices.
pub(crate) fn contiguous_index_span(indices: &[usize]) -> Option<(usize, usize)> {
    let first_index = *indices.first()?;
    for (offset, &index) in indices.iter().enumerate() {
        if index != first_index + offset {
            return None;
        }
    }
    Some((first_index, first_index + indices.len()))
}

#[cfg(test)]
mod tests {
    include!("common_tests.rs");
}
