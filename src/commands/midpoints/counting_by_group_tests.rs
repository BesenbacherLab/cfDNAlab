use super::*;
use crate::shared::length_axis::LengthAxis;
use ndarray::Array1;
use ndarray_npy::NpzWriter;
use std::fs::File;
use std::path::PathBuf;
use std::sync::Arc;

fn length_axis() -> Arc<LengthAxis> {
    Arc::new(LengthAxis::new(vec![20, 50, 100]).expect("test length axis should be valid"))
}

fn make_dense_counts() -> ProfileGroupsCounts {
    ProfileGroupsCounts::new(5, 3, length_axis())
}

fn make_sparse_counts() -> SparseProfileGroupsCounts {
    SparseProfileGroupsCounts::new(5, 3, length_axis())
}

fn write_sparse_partial_file(
    path: &Path,
    idx: &[u64],
    data: &[f32],
    shape: &[u64],
) -> Result<()> {
    let file = File::create(path)?;
    let mut npz = NpzWriter::new(file);
    npz.add_array("idx", &Array1::from(idx.to_vec()))?;
    npz.add_array("data", &Array1::from(data.to_vec()))?;
    npz.add_array("shape", &Array1::from(shape.to_vec()))?;
    npz.finish()?;
    Ok(())
}

fn assert_approx_eq(observed: f32, expected: f32) {
    assert!(
        (observed - expected).abs() <= 1e-6,
        "expected {expected}, got {observed}"
    );
}

#[test]
fn sparse_and_dense_indexing_match() -> Result<()> {
    let dense = make_dense_counts();
    let sparse = make_sparse_counts();

    for (position, group_idx, length) in [
        (0_usize, 0_usize, 20_usize),
        (3, 1, 49),
        (3, 1, 50),
        (4, 2, 99),
    ] {
        assert_eq!(
            dense.index_of(position, group_idx, length)?,
            sparse.index_of(position, group_idx, length)?
        );
    }

    assert!(dense.index_of(0, 0, 100).is_err());
    assert!(sparse.index_of(0, 0, 100).is_err());
    assert!(dense.index_of(0, 0, 19).is_err());
    assert!(sparse.index_of(0, 0, 19).is_err());
    assert!(dense.index_of(5, 0, 20).is_err());
    assert!(sparse.index_of(5, 0, 20).is_err());
    assert!(dense.index_of(0, 3, 20).is_err());
    assert!(sparse.index_of(0, 3, 20).is_err());

    assert_eq!(dense.index_of(0, 0, 20)?, 0);
    assert_eq!(dense.index_of(4, 0, 20)?, 4);
    assert_eq!(dense.index_of(0, 0, 50)?, 5);
    assert_eq!(dense.index_of(0, 1, 20)?, 10);

    Ok(())
}

#[test]
fn sparse_increment_accumulates_duplicate_entries() -> Result<()> {
    let mut sparse = make_sparse_counts();

    sparse.incr_weighted(2, 1, 49, 1.25)?;
    sparse.incr_weighted(2, 1, 49, 2.75)?;

    let flat_idx = sparse.index_of(2, 1, 49)?;
    assert_eq!(sparse.counts.len(), 1);
    assert_approx_eq(sparse.counts[&flat_idx], 4.0);

    Ok(())
}

#[test]
fn sparse_temp_write_read_roundtrip_keeps_sorted_indices_values_and_shape() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;
    let path = temp_dir.path().join("partial_file.npz");
    let mut sparse = make_sparse_counts();

    sparse.incr_weighted(4, 2, 99, 3.5)?;
    sparse.incr_weighted(0, 0, 20, 1.0)?;
    sparse.incr_weighted(3, 1, 50, 2.25)?;
    sparse.write_npz(&path)?;

    let dense = make_dense_counts();
    let partial_file = read_sparse_profile_partial_file(&path, [3, 2, 5], dense.counts.len())?;
    let expected_indices = vec![
        dense.index_of(0, 0, 20)? as u64,
        dense.index_of(3, 1, 50)? as u64,
        dense.index_of(4, 2, 99)? as u64,
    ];

    assert_eq!(partial_file.idx, expected_indices);
    assert_eq!(partial_file.data, vec![1.0, 2.25, 3.5]);

    Ok(())
}

#[test]
fn sparse_partial_file_reader_rejects_malformed_partial_files() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;
    let dense = make_dense_counts();
    let expected_shape = [3, 2, 5];
    let destination_len = dense.counts.len();

    let mismatched_lengths = temp_dir.path().join("mismatched_lengths.npz");
    write_sparse_partial_file(&mismatched_lengths, &[0, 1], &[1.0], &expected_shape)?;
    assert!(
        read_sparse_profile_partial_file(&mismatched_lengths, expected_shape, destination_len)
            .is_err()
    );

    let wrong_shape = temp_dir.path().join("wrong_shape.npz");
    write_sparse_partial_file(&wrong_shape, &[0], &[1.0], &[3, 3, 5])?;
    assert!(
        read_sparse_profile_partial_file(&wrong_shape, expected_shape, destination_len).is_err()
    );

    let descending = temp_dir.path().join("descending.npz");
    write_sparse_partial_file(&descending, &[2, 1], &[1.0, 1.0], &expected_shape)?;
    assert!(
        read_sparse_profile_partial_file(&descending, expected_shape, destination_len).is_err()
    );

    let out_of_bounds = temp_dir.path().join("out_of_bounds.npz");
    write_sparse_partial_file(
        &out_of_bounds,
        &[u64::try_from(destination_len)?],
        &[1.0],
        &expected_shape,
    )?;
    assert!(
        read_sparse_profile_partial_file(&out_of_bounds, expected_shape, destination_len).is_err()
    );

    Ok(())
}

#[test]
fn parallel_sparse_merge_sums_overlaps_across_chunks() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;
    let first_path = temp_dir.path().join("first.npz");
    let second_path = temp_dir.path().join("second.npz");
    let mut dense = ProfileGroupsCounts::new(
        4,
        2,
        Arc::new(LengthAxis::new(vec![20, 50, 100]).expect("test length axis should be valid")),
    );
    let shape = [2, 2, 4];
    let first_idx = dense.index_of(0, 0, 20)? as u64;
    let middle_idx = dense.index_of(1, 0, 50)? as u64;
    let last_idx = dense.index_of(3, 1, 50)? as u64;

    write_sparse_partial_file(&first_path, &[first_idx, last_idx], &[1.0, 2.0], &shape)?;
    write_sparse_partial_file(
        &second_path,
        &[first_idx, middle_idx, last_idx],
        &[3.0, 4.0, 5.0],
        &shape,
    )?;

    dense.add_from_sparse_npz_files_parallel_with_chunk_size(vec![first_path, second_path], 3)?;

    assert_approx_eq(dense.get(0, 0, 20)?, 4.0);
    assert_approx_eq(dense.get(1, 0, 50)?, 4.0);
    assert_approx_eq(dense.get(3, 1, 50)?, 7.0);
    assert_approx_eq(dense.get(0, 1, 20)?, 0.0);

    Ok(())
}

#[test]
fn parallel_sparse_merge_rejects_zero_chunk_size() {
    let mut dense = ProfileGroupsCounts::new(4, 2, length_axis());

    let error = dense
        .add_from_sparse_npz_files_parallel_with_chunk_size(Vec::<PathBuf>::new(), 0)
        .expect_err("zero chunk size should fail");

    assert!(
        error
            .to_string()
            .contains("sparse midpoint merge chunk size must be greater than zero"),
        "unexpected error: {error}"
    );
}

#[test]
fn sparse_merge_start_chunk_is_aligned_and_in_bounds() {
    let num_chunks = 30;
    let num_threads = rayon::current_num_threads().max(1);
    let start_stride = (num_chunks / (num_threads * 2)).max(1);

    for path_idx in 0..100 {
        let start_chunk = sparse_merge_start_chunk(path_idx, num_chunks);
        assert!(
            start_chunk < num_chunks,
            "start chunk {start_chunk} should be in bounds for {num_chunks} chunks"
        );
        assert_eq!(
            start_chunk % start_stride,
            0,
            "start chunk {start_chunk} should align to stride {start_stride}"
        );
    }

    assert_eq!(sparse_merge_start_chunk(10, 0), 0);
    assert_eq!(sparse_merge_start_chunk(10, 1), 0);
}
