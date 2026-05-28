use super::*;
use anyhow::Result;

fn codes_as_u64(codes: &KmerCodes) -> Vec<u64> {
    match codes {
        KmerCodes::U8(values) => values.iter().map(|&value| value as u64).collect(),
        KmerCodes::U16(values) => values.iter().map(|&value| value as u64).collect(),
        KmerCodes::U32(values) => values.iter().map(|&value| value as u64).collect(),
        KmerCodes::U64(values) => values.clone(),
    }
}

fn acgt_kmer_from_index(mut code_index: usize, k: usize) -> Vec<u8> {
    let mut kmer = vec![b'A'; k];
    for position in (0..k).rev() {
        kmer[position] = BASES[code_index % 4] as u8;
        code_index /= 4;
    }
    kmer
}

#[test]
fn subspace_spec_reuses_radix_encoding_when_k_fits() -> Result<()> {
    // Arrange: two selected 2-mers should be represented by compact selected codes 0 and 1.
    let selected = vec![&b"AC"[..], &b"GT"[..]];

    // Act
    let spec = build_subspace_kmer_spec(2, &selected)?;

    // Assert: small k reuses the existing radix-5 encoder and stores compact codes as u8.
    assert!(matches!(spec.encoding, SubspaceKmerEncoding::Radix5 { .. }));
    assert_eq!(spec.sentinel_missing(), u8::MAX as u64);
    assert_eq!(spec.encode_kmer_bytes(b"AC"), 0);
    assert_eq!(spec.encode_kmer_bytes(b"gt"), 1);
    assert_eq!(spec.encode_kmer_bytes(b"CG"), spec.sentinel_missing());
    assert_eq!(spec.encode_kmer_bytes(b"AN"), spec.sentinel_missing());

    // Derivation for ACGTAC with k=2:
    // AC -> 0, CG -> missing, GT -> 1, TA -> missing, AC -> 0, tail -> missing.
    let codes = spec.build_left_aligned_codes(b"ACGTAC");
    assert!(matches!(codes, KmerCodes::U8(_)));
    assert_eq!(
        codes_as_u64(&codes),
        vec![0, u8::MAX as u64, 1, u8::MAX as u64, 0, u8::MAX as u64]
    );

    Ok(())
}

#[test]
fn subspace_spec_uses_byte_lookup_when_k_exceeds_radix_limit() -> Result<()> {
    // Arrange: 28-mers do not fit the existing radix-5 u64 representation.
    let first = b"ACGT".repeat(7);
    let lowercase_first: Vec<u8> = first
        .iter()
        .map(|base| base.to_ascii_lowercase())
        .collect();
    let second = b"TGCA".repeat(7);
    let selected = vec![first.as_slice(), lowercase_first.as_slice(), second.as_slice()];

    // Act
    let spec = build_subspace_kmer_spec(28, &selected)?;

    // Assert: the selected-code contract is the same even though the backing lookup is byte-based,
    // and the duplicate lowercase first k-mer does not shift the second unique code.
    assert!(matches!(spec.encoding, SubspaceKmerEncoding::Bytes { .. }));
    assert_eq!(spec.sentinel_missing(), u8::MAX as u64);
    assert_eq!(spec.encode_kmer_bytes(first.as_slice()), 0);
    assert_eq!(spec.encode_kmer_bytes(second.as_slice()), 1);
    assert_eq!(spec.encode_kmer_bytes(&lowercase_first), 0);

    // Derivation: first 28 bases match selected code 0. The shifted 28-mer contains 27 bases from
    // `first` plus T and is not selected. The remaining 27 positions cannot start a full k-mer.
    let mut reference = first.clone();
    reference.push(b'T');
    let codes = spec.build_left_aligned_codes(&reference);
    let mut expected = vec![spec.sentinel_missing(); reference.len()];
    expected[0] = 0;
    assert!(matches!(codes, KmerCodes::U8(_)));
    assert_eq!(codes_as_u64(&codes), expected);

    Ok(())
}

#[test]
fn subspace_width_uses_smallest_dtype_with_missing_sentinel() -> Result<()> {
    // The selected code range is 0..n_selected, and the top value of the dtype is reserved for
    // missing. Therefore 254 selected kmers still fit in u8, while 255 requires u16.
    assert_eq!(choose_subspace_width(254)?, (Width::U8, u8::MAX as u64));
    assert_eq!(
        choose_subspace_width(255)?,
        (Width::U16, u16::MAX as u64)
    );
    assert_eq!(
        choose_subspace_width(65_535)?,
        (Width::U32, u32::MAX as u64)
    );

    Ok(())
}

#[test]
fn subspace_spec_sizes_storage_from_unique_normalized_kmers() -> Result<()> {
    // Arrange: 254 unique 4-mers fit in u8 when the missing sentinel uses 255. Adding a duplicate
    // raw entry brings the input length to 255, which would require u16 if sizing happened before
    // deduplication.
    let mut selected: Vec<Vec<u8>> = (0..254)
        .map(|code_index| acgt_kmer_from_index(code_index, 4))
        .collect();
    selected.push(b"aaaa".to_vec());

    // Act
    let spec = build_subspace_kmer_spec(4, &selected)?;

    // Assert
    assert_eq!(spec.sentinel_missing(), u8::MAX as u64);
    assert_eq!(spec.encode_kmer_bytes(b"AAAA"), 0);
    assert_eq!(spec.encode_kmer_bytes(b"aaaa"), 0);
    assert_eq!(spec.encode_kmer_bytes(selected[253].as_slice()), 253);

    Ok(())
}

#[test]
fn subspace_spec_dedupes_kmers_after_normalization() -> Result<()> {
    // Arrange: AC and ac normalize to the same selected k-mer. GT should still receive code 1,
    // proving that duplicate normalized entries do not consume selected-code space.
    let selected = vec![&b"AC"[..], &b"ac"[..], &b"GT"[..]];

    // Act
    let spec = build_subspace_kmer_spec(2, &selected)?;

    // Assert
    assert_eq!(spec.encode_kmer_bytes(b"AC"), 0);
    assert_eq!(spec.encode_kmer_bytes(b"ac"), 0);
    assert_eq!(spec.encode_kmer_bytes(b"GT"), 1);

    Ok(())
}

#[test]
fn subspace_spec_rejects_invalid_selected_kmers() {
    // Arrange + Act
    let selected = vec![&b"AN"[..]];
    let error =
        build_subspace_kmer_spec(2, &selected).expect_err("N-containing selected k-mer should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid selected k-mer"),
        "unexpected error: {error}"
    );
}
