use super::FragmentLengthArgs;

#[test]
fn fragment_length_args_rejects_inverted_range() {
    let args = FragmentLengthArgs {
        min_fragment_length: 500,
        max_fragment_length: 100,
    };

    let error = args
        .validate()
        .expect_err("inverted fragment length range should fail");
    let message = error.to_string();

    assert!(
        message.contains("--min-fragment-length (500) must be <= --max-fragment-length (100)"),
        "unexpected error: {message}"
    );
}

#[test]
fn fragment_length_args_accepts_inclusive_single_length_range() {
    let args = FragmentLengthArgs {
        min_fragment_length: 10,
        max_fragment_length: 10,
    };

    args.validate()
        .expect("single-length fragment range should be valid");
}
