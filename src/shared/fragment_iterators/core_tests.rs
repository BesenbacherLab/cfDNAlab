use super::ensure_nondecreasing_bam_coordinates;
use rust_htslib::bam::Record;

fn record_at(tid: i32, pos: i64) -> Record {
    let mut rec = Record::new();
    rec.set_tid(tid);
    rec.set_pos(pos);
    rec
}

#[test]
fn ensure_nondecreasing_bam_coordinates_accepts_equal_or_increasing_positions() {
    let mut last = None;

    ensure_nondecreasing_bam_coordinates(&mut last, &record_at(0, 10))
        .expect("first record should be accepted");
    ensure_nondecreasing_bam_coordinates(&mut last, &record_at(0, 10))
        .expect("equal position should be accepted");
    ensure_nondecreasing_bam_coordinates(&mut last, &record_at(0, 11))
        .expect("larger position should be accepted");
}

#[test]
fn ensure_nondecreasing_bam_coordinates_rejects_backwards_positions_on_same_tid() {
    let mut last = Some((0, 20));

    let error = ensure_nondecreasing_bam_coordinates(&mut last, &record_at(0, 19))
        .expect_err("backwards read position should fail");
    assert!(
        error.to_string().contains(
            "coordinate-sorted with nondecreasing read.pos within single-chromosome stream",
        ),
        "unexpected error message: {error}"
    );
}

#[test]
fn ensure_nondecreasing_bam_coordinates_rejects_tid_changes_within_one_stream() {
    let mut last = Some((1, 5));

    let error = ensure_nondecreasing_bam_coordinates(&mut last, &record_at(0, 100))
        .expect_err("tid change should fail");
    assert!(
        error
            .to_string()
            .contains("multiple tids inside single-chromosome stream"),
        "unexpected error message: {error}"
    );
}
