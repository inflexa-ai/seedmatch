//! Search rule 2: the windows of the pattern length with `min_mm` to `max_mm` mismatches.
//!
//! The 1-based starts go from 2 - m to n. Each pattern letter outside the sequence is a mismatch.

use std::iter;

use crate::letters::{Letters, WORD};
use crate::{to_i64, to_u64};

/// Find the matches of rule 2 in `seq`. `emit` gets the 1-based start of each match, in
/// increasing order.
pub fn search(seq: &[u8], letters: &Letters, max_mm: u64, min_mm: u64, emit: &mut impl FnMut(i64)) {
    if letters.len() <= WORD {
        bits(seq, letters, max_mm, min_mm, emit);
    } else {
        scalar(seq, letters, max_mm, min_mm, emit);
    }
}

/// Shift-and with k mismatches, for a pattern of 1 to `WORD` letters.
///
/// Bit a of state d is 1 if the pattern letters 0 to a match the last a + 1 letters with d
/// mismatches or fewer. The text has m - 1 virtual letters before and after the sequence, and a
/// virtual letter matches no pattern letter. Thus the windows of the text are the windows of
/// rule 2, and the result is exact.
pub fn bits(seq: &[u8], letters: &Letters, max_mm: u64, min_mm: u64, emit: &mut impl FnMut(i64)) {
    let m = letters.len();
    debug_assert!((1..=WORD).contains(&m));
    let top = m - 1;
    // A window has m mismatches or fewer. Thus if max_mm >= m, each window has max_mm mismatches
    // or fewer, and no state is necessary for this limit.
    let upper = usize::try_from(max_mm).ok().filter(|&max_mm| max_mm < m);
    // The window has min_mm mismatches or more if state min_mm - 1 does not match it.
    let below = match usize::try_from(min_mm) {
        Ok(0) => None,
        Ok(min_mm) if min_mm <= m => Some(min_mm - 1),
        _ => return,
    };
    let Some(last) = upper.max(below) else {
        for start in 2 - to_i64(m)..=to_i64(seq.len()) {
            emit(start);
        }
        return;
    };
    let words = letters.words();
    let mut states = [0_u64; WORD];
    let states = &mut states[..=last];
    // After the m - 1 virtual letters before the sequence, state d matches each prefix of d
    // letters or fewer.
    for (d, state) in states.iter_mut().enumerate() {
        *state = (1 << d.min(top)) - 1;
    }
    let text = seq.iter().map(|&code| words[usize::from(code & 15)]);
    // The window that ends at virtual letter v after the sequence has v mismatches or more. Thus
    // only the first max_mm of these letters can end a match.
    let after = iter::repeat_n(0, upper.map_or(top, |max_mm| max_mm.min(top)));
    let ends = |state: u64| (state >> top) & 1 == 1;
    for (start, word) in (2 - to_i64(m)..).zip(text.chain(after)) {
        step(states, word);
        if upper.is_none_or(|upper| ends(states[upper]))
            && below.is_none_or(|below| !ends(states[below]))
        {
            emit(start);
        }
    }
}

/// Read one text letter, with the word of its mask.
fn step(states: &mut [u64], word: u64) {
    let mut before = states[0];
    states[0] = ((before << 1) | 1) & word;
    for state in &mut states[1..] {
        let old = *state;
        *state = (((old << 1) | 1) & word) | ((before << 1) | 1);
        before = old;
    }
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
