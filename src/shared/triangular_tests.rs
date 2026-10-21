use super::*;

#[test]
fn triangular_weights_rise_to_the_center_and_fall_symmetrically() {
    let cache = TriangularWeightCache::default();
    for (radius, expected) in [
        (0, vec![1]),
        (1, vec![1, 2, 1]),
        (2, vec![1, 2, 3, 2, 1]),
        (7, vec![1, 2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 1]),
    ] {
        let weights = cache.get(radius).unwrap();

        assert_eq!(weights.as_ref(), expected.as_slice());
    }
}

#[test]
fn clipped_triangles_retain_the_original_peak_position_and_slope() {
    let cache = TriangularWeightCache::default();
    // A radius-2 triangle is [1,2,3,2,1]. At the left end skip its unavailable
    // two neighbors, and at the right end keep only the first three weights
    let weights = cache.get(2).unwrap();
    // A radius-3 triangle is [1,2,3,4,3,2,1]. One left neighbor gives [3,4,3,2,1]
    let asymmetric = cache.get(3).unwrap();

    assert_eq!(&weights[2..], [3, 2, 1]);
    assert_eq!(&weights[..3], [1, 2, 3]);
    assert_eq!(&asymmetric[2..], [3, 4, 3, 2, 1]);
}

#[test]
fn different_radii_are_cached_independently_and_reused_when_requested_again() {
    let cache = TriangularWeightCache::default();

    let larger = cache.get(2).unwrap();
    let smaller = cache.get(1).unwrap();
    let larger_again = cache.get(2).unwrap();
    let smaller_again = cache.get(1).unwrap();

    assert!(Arc::ptr_eq(&larger, &larger_again));
    assert!(Arc::ptr_eq(&smaller, &smaller_again));
    assert!(!Arc::ptr_eq(&larger, &smaller));
    // A separate run owns separate profiles, rather than extending a process-wide cache
    let separate_run = TriangularWeightCache::default();
    assert!(!Arc::ptr_eq(&larger, &separate_run.get(2).unwrap()));
}

#[test]
fn parallel_requests_share_the_same_cached_profile() {
    let cache = TriangularWeightCache::default();
    let start = std::sync::Barrier::new(4);

    let profiles = std::thread::scope(|scope| {
        let handles = (0..4)
            .map(|_| {
                scope.spawn(|| {
                    start.wait();
                    cache.get(2).unwrap()
                })
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .map(|handle| handle.join().unwrap())
            .collect::<Vec<_>>()
    });

    assert_eq!(profiles.len(), 4);
    for profile in &profiles[1..] {
        assert!(Arc::ptr_eq(&profiles[0], profile));
    }
}

#[test]
fn invalid_kernel_sizes_return_errors_without_poisoning_the_cache() {
    let cache = TriangularWeightCache::default();
    for (radius, expected_message) in [
        (usize::MAX, "peak exceeds usize"),
        (usize::MAX / 2 + 1, "length exceeds usize"),
        // MAX elements cannot fit in Vec<usize>'s byte capacity, so reservation fails
        // before any allocation or iteration over this deliberately impossible profile
        (usize::MAX / 2, "allocating triangular weights"),
    ] {
        let error = cache.get(radius).unwrap_err();
        assert!(error.to_string().contains(expected_message));
    }

    assert_eq!(cache.get(1).unwrap().as_ref(), [1, 2, 1]);
}
