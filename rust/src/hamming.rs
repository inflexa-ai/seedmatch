//! Search rule 2: the windows of the pattern length with `min_mm` to `max_mm` mismatches.
//!
//! The 1-based starts go from 2 - m to n. Each pattern letter outside the sequence is a mismatch.

use crate::letters::Letters;
use crate::{to_i64, to_u64};

/// Find the matches of rule 2 in `seq`. `emit` gets the 1-based start of each match, in
/// increasing order.
pub fn search(seq: &[u8], letters: &Letters, max_mm: u64, min_mm: u64, emit: &mut impl FnMut(i64)) {
    scalar(seq, letters, max_mm, min_mm, emit);
}

/// Count the mismatches of each window, with an early stop.
pub fn scalar(seq: &[u8], letters: &Letters, max_mm: u64, min_mm: u64, emit: &mut impl FnMut(i64)) {
    let is_match = |count: u64| min_mm <= count && count <= max_mm;
    for lead in (1..letters.len()).rev() {
        if is_match(mismatches(seq, letters, lead, 0, max_mm)) {
            emit(1 - to_i64(lead));
        }
    }
    for at in 0..seq.len() {
        if is_match(mismatches(seq, letters, 0, at, max_mm)) {
            emit(to_i64(at) + 1);
        }
    }
}

/// Count the mismatches of the window where pattern letter `first` is on `seq[at]`, and the
/// letters before `first` are before the sequence. The count stops after it passes `max_mm`.
fn mismatches(seq: &[u8], letters: &Letters, first: usize, at: usize, max_mm: u64) -> u64 {
    let m = letters.len();
    let inside = (m - first).min(seq.len() - at);
    let mut count = to_u64(m - inside);
    for (a, &code) in (first..).zip(&seq[at..at + inside]) {
        if count > max_mm {
            break;
        }
        if !letters.ok(a, code) {
            count += 1;
        }
    }
    count
}
