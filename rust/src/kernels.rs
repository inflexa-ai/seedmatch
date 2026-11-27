//! The kernels on a set of sequences: the arrays `codes` and `offsets` of a `SequenceSet`.
//!
//! Sequence i is `codes[offsets[i]..offsets[i + 1]]`. The kernels give the sequences in the order
//! of the set, and the matches of each sequence in the order of Biostrings. The coordinates are
//! 1-based and inclusive.

use std::fmt;

use crate::letters::{Compat, Letters};
use crate::{hamming, indel, to_i64};

/// The sequence index, the 1-based start, and the 1-based end of each match.
pub type Ranges = (Vec<i64>, Vec<i64>, Vec<i64>);

/// The offsets do not give sequences inside `codes`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OffsetsError {
    pub codes: usize,
}

impl fmt::Display for OffsetsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "offsets do not fit codes of length {}; expected the offsets of a SequenceSet",
            self.codes
        )
    }
}

impl std::error::Error for OffsetsError {}

/// Call `visit` with the index and the codes of each sequence.
///
/// A write through a base array of the `SequenceSet` can change the offsets during the search.
/// Thus each offset is read one time, and checked before it is used.
fn each_sequence(
    codes: &[u8],
    offsets: &[i64],
    mut visit: impl FnMut(usize, &[u8]),
) -> Result<(), OffsetsError> {
    let error = OffsetsError { codes: codes.len() };
    let Some((&first, rest)) = offsets.split_first() else {
        return Ok(());
    };
    let mut lo = usize::try_from(first).map_err(|_| error)?;
    if lo > codes.len() {
        return Err(error);
    }
    for (i, &offset) in rest.iter().enumerate() {
        let hi = usize::try_from(offset).map_err(|_| error)?;
        if hi < lo || hi > codes.len() {
            return Err(error);
        }
        visit(i, &codes[lo..hi]);
        lo = hi;
    }
    Ok(())
}

/// Return the max_mm of the indel DP, or `None` if the search gives the result of rule 2: rule 5
/// for max_mm 0, and rule 4b for max_mm at or above the pattern length.
fn indel_budget(max_mm: u64, letters: &Letters) -> Option<usize> {
    let budget = usize::try_from(max_mm).ok()?;
    (budget > 0 && budget < letters.len()).then_some(budget)
}

fn counts_of(
    codes: &[u8],
    offsets: &[i64],
    mut search: impl FnMut(&[u8]) -> i64,
) -> Result<Vec<i64>, OffsetsError> {
    let mut counts = Vec::with_capacity(offsets.len().saturating_sub(1));
    each_sequence(codes, offsets, |_, seq| counts.push(search(seq)))?;
    Ok(counts)
}

fn hamming_counts_of(
    codes: &[u8],
    offsets: &[i64],
    letters: &Letters,
    max_mm: u64,
    min_mm: u64,
) -> Result<Vec<i64>, OffsetsError> {
    counts_of(codes, offsets, |seq| {
        let mut found = 0;
        hamming::search(seq, letters, max_mm, min_mm, &mut |_| found += 1);
        found
    })
}

fn hamming_ranges_of(
    codes: &[u8],
    offsets: &[i64],
    letters: &Letters,
    max_mm: u64,
    min_mm: u64,
) -> Result<(Vec<i64>, Vec<i64>), OffsetsError> {
    let mut index = Vec::new();
    let mut starts = Vec::new();
    each_sequence(codes, offsets, |i, seq| {
        let i = to_i64(i);
        hamming::search(seq, letters, max_mm, min_mm, &mut |start| {
            index.push(i);
            starts.push(start);
        });
    })?;
    Ok((index, starts))
}

/// Return the number of matches of rule 2 in each sequence.
pub fn hamming_counts(
    codes: &[u8],
    offsets: &[i64],
    pattern: &[u8],
    compat: &Compat,
    max_mm: u64,
    min_mm: u64,
) -> Result<Vec<i64>, OffsetsError> {
    hamming_counts_of(
        codes,
        offsets,
        &Letters::new(pattern, compat),
        max_mm,
        min_mm,
    )
}

/// Return the sequence index and the 1-based start of each match of rule 2. The end of a match is
/// its start + len(pattern) - 1.
pub fn hamming_ranges(
    codes: &[u8],
    offsets: &[i64],
    pattern: &[u8],
    compat: &Compat,
    max_mm: u64,
    min_mm: u64,
) -> Result<(Vec<i64>, Vec<i64>), OffsetsError> {
    hamming_ranges_of(
        codes,
        offsets,
        &Letters::new(pattern, compat),
        max_mm,
        min_mm,
    )
}

/// Return the number of matches of the indel search in each sequence (rules 4, 4b, and 5).
pub fn indel_counts(
    codes: &[u8],
    offsets: &[i64],
    pattern: &[u8],
    compat: &Compat,
    max_mm: u64,
) -> Result<Vec<i64>, OffsetsError> {
    let letters = Letters::new(pattern, compat);
    let Some(budget) = indel_budget(max_mm, &letters) else {
        return hamming_counts_of(codes, offsets, &letters, max_mm, 0);
    };
    counts_of(codes, offsets, |seq| {
        let mut found = 0;
        indel::search(seq, &letters, budget, &mut |_, _| found += 1);
        found
    })
}

/// Return the sequence index, the 1-based start, and the 1-based end of each match of the indel
/// search (rules 4, 4b, and 5).
pub fn indel_ranges(
    codes: &[u8],
    offsets: &[i64],
    pattern: &[u8],
    compat: &Compat,
    max_mm: u64,
) -> Result<Ranges, OffsetsError> {
    let letters = Letters::new(pattern, compat);
    let Some(budget) = indel_budget(max_mm, &letters) else {
        let (index, starts) = hamming_ranges_of(codes, offsets, &letters, max_mm, 0)?;
        let last = to_i64(letters.len()) - 1;
        let ends = starts.iter().map(|start| start + last).collect();
        return Ok((index, starts, ends));
    };
    let mut index = Vec::new();
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    each_sequence(codes, offsets, |i, seq| {
        let i = to_i64(i);
        indel::search(seq, &letters, budget, &mut |start, end| {
            index.push(i);
            starts.push(start);
            ends.push(end);
        });
    })?;
    Ok((index, starts, ends))
}
