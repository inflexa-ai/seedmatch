//! Search rule 4: the matches with indels, and one held candidate.
//!
//! An indel match stays inside its sequence. The caller makes sure that 1 <= max_mm < len().

use crate::letters::{Letters, WORD};
use crate::to_i64;

/// The positions of one block of the filter. The candidates of a block take `BLOCK` bits, thus the
/// memory of the filter does not grow with the sequence.
pub const BLOCK: usize = 1 << 16;

/// Find the matches of rule 4 in `seq`. `emit` gets the 1-based start and end of each match.
pub fn search(seq: &[u8], letters: &Letters, max_mm: usize, emit: &mut impl FnMut(i64, i64)) {
    if letters.len() <= WORD {
        bits(seq, letters, max_mm, BLOCK, emit);
    } else {
        scalar(seq, letters, max_mm, emit);
    }
}

/// Run the exact step at each position of `seq`.
pub fn scalar(seq: &[u8], letters: &Letters, max_mm: usize, emit: &mut impl FnMut(i64, i64)) {
    let mut row = dp_row(letters, max_mm);
    let mut held = Held::default();
    for j in 0..seq.len() {
        if let Some((end, cost)) = accept(seq, letters, max_mm, j, &mut row) {
            held.offer(j + 1, end, cost, emit);
        }
    }
    held.finish(emit);
}

/// Run the exact step only at the candidates of the filter, for a pattern of 1 to `WORD` letters.
///
/// Each position that the exact step accepts is a candidate. At a position that it does not
/// accept, `scalar` changes nothing. Thus the result is the result of `scalar`.
pub fn bits(
    seq: &[u8],
    letters: &Letters,
    max_mm: usize,
    block: usize,
    emit: &mut impl FnMut(i64, i64),
) {
    let mut row = dp_row(letters, max_mm);
    let mut held = Held::default();
    filter(seq, letters, max_mm, block, |j| {
        if let Some((end, cost)) = accept(seq, letters, max_mm, j, &mut row) {
            held.offer(j + 1, end, cost, emit);
        }
    });
    held.finish(emit);
}

/// Return each 0-based position that the exact step accepts, in increasing order.
pub fn accepted(seq: &[u8], letters: &Letters, max_mm: usize) -> Vec<usize> {
    let mut row = dp_row(letters, max_mm);
    (0..seq.len())
        .filter(|&j| accept(seq, letters, max_mm, j, &mut row).is_some())
        .collect()
}

/// Return each 0-based candidate of the filter, in increasing order.
pub fn candidates(seq: &[u8], letters: &Letters, max_mm: usize, block: usize) -> Vec<usize> {
    let mut found = Vec::new();
    filter(seq, letters, max_mm, block, |j| found.push(j));
    found
}

/// Call `visit` for each position j with S(j) <= max_mm, in increasing order. S(j) is the least
/// edit distance of the pattern and a text part that starts at j and ends inside `seq`.
///
/// The Myers (1999) bit-vector algorithm gives, for each end of a text part, the least edit
/// distance over each start. On the reversed pattern and the reversed text, each step thus gives
/// S(j) for one j.
fn filter(
    seq: &[u8],
    letters: &Letters,
    max_mm: usize,
    block: usize,
    mut visit: impl FnMut(usize),
) {
    let m = letters.len();
    debug_assert!((1..=WORD).contains(&m) && block > 0);
    let peq = letters.reversed_words();
    let high = 1_u64 << (m - 1);
    let n = seq.len();
    let mut marks = vec![0_u64; block.min(n).div_ceil(WORD)];
    let mut b0 = 0;
    while b0 < n {
        let b1 = n.min(b0 + block);
        // A text part with an edit distance of max_mm or less has m + max_mm letters or fewer.
        // Thus each such part that starts in the block ends before `top`, and the run from `top`
        // gives the exact decision S(j) <= max_mm for each j of the block.
        let top = n.min(b1 + m + max_mm - 1);
        marks.fill(0);
        let (mut pv, mut mv, mut score) = (u64::MAX, 0_u64, m);
        for i in (b0..top).rev() {
            let eq = peq[usize::from(seq[i] & 15)];
            let xv = eq | mv;
            let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
            let mut ph = mv | !(xh | pv);
            let mut mh = pv & xh;
            if ph & high != 0 {
                score += 1;
            }
            if mh & high != 0 {
                score -= 1;
            }
            // No `| 1`: a text part can start at each position of the run.
            ph <<= 1;
            mh <<= 1;
            pv = mh | !(xv | ph);
            mv = ph & xv;
            if i < b1 && score <= max_mm {
                let k = i - b0;
                marks[k / WORD] |= 1 << (k % WORD);
            }
        }
        for (w, &word) in marks.iter().enumerate() {
            let mut rest = word;
            while rest != 0 {
                visit(b0 + w * WORD + rest.trailing_zeros() as usize);
                rest &= rest - 1;
            }
        }
        b0 = b1;
    }
}

/// A candidate aligns len() - 1 letters or fewer with a budget of max_mm or less. Thus its widths
/// are 0 to len() - 1 + max_mm.
fn dp_row(letters: &Letters, max_mm: usize) -> Vec<usize> {
    vec![0; letters.len() + max_mm]
}

/// Return the 1-based end and the cost of the match that starts at the 0-based position `j`, or
/// `None` if the cost is more than `max_mm`.
///
/// The pattern letters before the first letter that matches `seq[j]` are deletions. The DP aligns
/// the letters after it to the text after `j`.
fn accept(
    seq: &[u8],
    letters: &Letters,
    max_mm: usize,
    j: usize,
    row: &mut [usize],
) -> Option<(usize, usize)> {
    let before = letters.first(seq[j]);
    if before > max_mm {
        return None;
    }
    let (distance, width) = right_alignment(seq, j + 1, letters, before + 1, max_mm - before, row)?;
    Some((j + 1 + width, before + distance))
}

/// Align the pattern letters from `letter` to `seq[begin..begin + w]`, for each width w.
///
/// Return the least edit distance and the shortest width with it, or `None` if the least distance
/// is more than `budget`.
fn right_alignment(
    seq: &[u8],
    begin: usize,
    letters: &Letters,
    letter: usize,
    budget: usize,
    row: &mut [usize],
) -> Option<(usize, usize)> {
    let rest = letters.len() - letter;
    let widths = (seq.len() - begin).min(rest + budget);
    let cap = budget + 1;
    // After k letters, the distance at width w is |k - w| or more. Thus only the band of widths
    // k - budget to k + budget can have a distance of budget or less, and the DP computes only that
    // band. Each value saturates at cap. A value of budget or less comes only from values of budget
    // or less, so the saturation keeps it exact. Each value of the row is written before it is
    // read, thus the result does not depend on the positions that the caller visited before.
    let mut low = 0;
    let high = widths.min(budget);
    for (w, value) in row[..=high].iter_mut().enumerate() {
        *value = w;
    }
    let mut best = 0;
    for k in 1..=rest {
        low = k.saturating_sub(budget);
        let high = widths.min(k + budget);
        if low > high {
            return None;
        }
        if high == k + budget {
            // Row k - 1 did not compute this width.
            row[high] = cap;
        }
        let mut diagonal;
        let mut left;
        if low == 0 {
            diagonal = row[0];
            row[0] = k;
            left = k;
            best = k;
        } else {
            diagonal = row[low - 1];
            left = cap;
            best = cap;
        }
        for w in low.max(1)..=high {
            let up = row[w];
            let substitution = if letters.ok(letter + k - 1, seq[begin + w - 1]) {
                diagonal
            } else {
                diagonal + 1
            };
            let value = substitution.min(up + 1).min(left + 1).min(cap);
            diagonal = up;
            row[w] = value;
            left = value;
            best = best.min(value);
        }
        if best > budget {
            return None;
        }
    }
    let width = low + row[low..].iter().position(|&value| value == best)?;
    Some((best, width))
}

struct Candidate {
    start: usize,
    end: usize,
    cost: usize,
}

/// The candidate that waits for the next accepted position. A later position with an end at or
/// before the held end replaces the held candidate if its cost is not more.
#[derive(Default)]
struct Held(Option<Candidate>);

impl Held {
    fn offer(&mut self, start: usize, end: usize, cost: usize, emit: &mut impl FnMut(i64, i64)) {
        if let Some(held) = &self.0
            && end <= held.end
        {
            if cost <= held.cost {
                self.0 = Some(Candidate { start, end, cost });
            }
            return;
        }
        if let Some(held) = self.0.replace(Candidate { start, end, cost }) {
            emit(to_i64(held.start), to_i64(held.end));
        }
    }

    fn finish(self, emit: &mut impl FnMut(i64, i64)) {
        if let Some(held) = self.0 {
            emit(to_i64(held.start), to_i64(held.end));
        }
    }
}
