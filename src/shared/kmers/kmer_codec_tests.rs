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
fn radix5_width_selection_uses_the_smallest_type_with_two_sentinels() -> Result<()> {
    // The full radix-5 code space needs two reserved sentinel values in the chosen integer type.
    // k = 3 fits in u8, k = 5 needs u16, k = 10 needs u32, and k = 26 needs u64.
    let (width_3, none_3, n_3) = choose_width(3)?;
    assert_eq!(width_3, Width::U8);
    assert_eq!(none_3, u8::MAX as u64);
    assert_eq!(n_3, (u8::MAX - 1) as u64);

    assert_eq!(choose_width(5)?.0, Width::U16);
    assert_eq!(choose_width(10)?.0, Width::U32);
    assert_eq!(choose_width(26)?.0, Width::U64);

    Ok(())
}

#[test]
fn radix5_left_aligned_codes_decode_each_position() -> Result<()> {
    // For ACGTAC and k = 2, starts 0..4 have complete 2-mers and the final position has no full
    // k-mer. The no-full-k-mer sentinel decodes to NN.
    let spec = build_kmer_specs(&[2])?
        .remove(&2)
        .expect("requested k-mer spec should exist");

    let codes = spec.build_left_aligned_codes(b"ACGTAC");
    let decoded: Vec<String> = codes.iter().map(|&code| spec.decode_kmer(code)).collect();

    assert_eq!(decoded, vec!["AC", "CG", "GT", "TA", "AC", "NN"]);
    Ok(())
}

#[test]
fn radix5_left_aligned_codes_mark_n_windows_and_tail_sentinels() -> Result<()> {
    // For ACGTACN and k = 3:
    //   starts 0..3 are ordinary ACG, CGT, GTA, TAC
    //   start 4 is ACN and must use the contains-N sentinel
    //   starts 5..6 have no full 3-mer and must use the no-full-k-mer sentinel.
    let spec = build_kmer_specs(&[3])?
        .remove(&3)
        .expect("requested k-mer spec should exist");

    let codes = spec.build_left_aligned_codes(b"ACGTACN");

    assert_eq!(spec.decode_kmer(codes[0]), "ACG");
    assert_eq!(spec.decode_kmer(codes[1]), "CGT");
    assert_eq!(spec.decode_kmer(codes[2]), "GTA");
    assert_eq!(spec.decode_kmer(codes[3]), "TAC");
    assert_eq!(codes[4], spec.sentinel_n());
    assert_eq!(codes[5], spec.sentinel_none());
    assert_eq!(codes[6], spec.sentinel_none());
    assert_eq!(spec.decode_kmer(codes[4]), "NNN");
    assert_eq!(spec.decode_kmer(codes[5]), "NNN");
    assert_eq!(spec.decode_kmer(codes[6]), "NNN");
    Ok(())
}

#[test]
fn subspace_spec_assigns_compact_selected_codes() -> Result<()> {
    // Arrange: two selected 2-mers should be represented by compact selected codes 0 and 1.
    let selected = vec![&b"AC"[..], &b"GT"[..]];

    // Act
    let spec = build_subspace_kmer_spec(2, &selected)?;

    // Assert: selected-subspace codes are compact and stored as u8.
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
fn subspace_spec_supports_k_above_the_radix_limit() -> Result<()> {
    // Arrange: 28-mers do not fit the existing radix-5 u64 representation.
    let first = b"ACGT".repeat(7);
    let lowercase_first: Vec<u8> = first.iter().map(|base| base.to_ascii_lowercase()).collect();
    let second = b"TGCA".repeat(7);
    let selected = vec![
        first.as_slice(),
        lowercase_first.as_slice(),
        second.as_slice(),
    ];

    // Act
    let spec = build_subspace_kmer_spec(28, &selected)?;

    // Assert: the selected-code contract still holds above the radix-5 limit, and the duplicate
    // lowercase first k-mer does not shift the second unique code.
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

    let lowercase_reference: Vec<u8> = reference
        .iter()
        .map(|base| base.to_ascii_lowercase())
        .collect();
    let lowercase_codes = spec.build_left_aligned_codes(&lowercase_reference);
    assert_eq!(codes_as_u64(&lowercase_codes), expected);

    Ok(())
}

#[test]
fn subspace_width_uses_smallest_dtype_with_missing_sentinel() -> Result<()> {
    // The selected code range is 0..n_selected, and the top value of the dtype is reserved for
    // missing. Therefore 254 selected kmers still fit in u8, while 255 requires u16.
    assert_eq!(choose_subspace_width(254)?, (Width::U8, u8::MAX as u64));
    assert_eq!(choose_subspace_width(255)?, (Width::U16, u16::MAX as u64));
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
    let error = build_subspace_kmer_spec(2, &selected)
        .expect_err("N-containing selected k-mer should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid selected k-mer"),
        "unexpected error: {error}"
    );
}
