//! The random input of the Rust tests.

use _kernels::letters::{Compat, MASKS};

/// The SplitMix64 generator: a fixed seed gives the same cases on each machine.
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Return a value from 0 to `bound` - 1.
    pub fn below(&mut self, bound: usize) -> usize {
        // The value is less than `bound`, thus the cast back to usize keeps it.
        (self.next_u64() % bound as u64) as usize
    }

    pub fn one_in(&mut self, odds: usize) -> bool {
        self.below(odds) == 0
    }
}

/// The table of `compat_table` in `src/seedmatch/search/compat.py`.
pub fn compat(fixed_pattern: bool, fixed_subject: bool) -> Compat {
    let mut table = [[false; MASKS]; MASKS];
    for (p, row) in table.iter_mut().enumerate().skip(1) {
        for (s, ok) in row.iter_mut().enumerate().skip(1) {
            *ok = match (fixed_pattern, fixed_subject) {
                (true, true) => p == s,
                (true, false) => p & !s == 0,
                (false, true) => s & !p == 0,
                (false, false) => p & s != 0,
            };
        }
    }
    table
}

/// Return a mask from 1 to 15. Three masks in four are a base: 1, 2, 4, or 8.
pub fn letter(rng: &mut SplitMix64) -> u8 {
    if rng.one_in(4) {
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15][rng.below(15)]
    } else {
        [1, 2, 4, 8][rng.below(4)]
    }
}

/// Return a sequence of `len` masks. One sequence in 10 is a run of N (15) or a run of A (1),
/// where most positions are candidates of the indel filter.
pub fn sequence(rng: &mut SplitMix64, len: usize) -> Vec<u8> {
    if rng.one_in(10) {
        let run = if rng.one_in(2) { 15 } else { 1 };
        return vec![run; len];
    }
    (0..len).map(|_| letter(rng)).collect()
}

/// Put a copy of the pattern with `edits` random edits at a random position of `seq`. Each letter
/// of the copy is a base of its pattern mask.
pub fn plant(rng: &mut SplitMix64, seq: &mut Vec<u8>, pattern: &[u8], edits: usize) {
    let mut copy: Vec<u8> = pattern.iter().map(|&mask| base_of(rng, mask)).collect();
    for _ in 0..edits {
        let at = rng.below(copy.len() + 1);
        match rng.below(3) {
            0 if at < copy.len() => copy[at] = letter(rng),
            1 => copy.insert(at, letter(rng)),
            _ if at < copy.len() => {
                copy.remove(at);
            }
            _ => copy.push(letter(rng)),
        }
    }
    let at = rng.below(seq.len() + 1);
    seq.splice(at..at, copy);
}

fn base_of(rng: &mut SplitMix64, mask: u8) -> u8 {
    let bases: Vec<u8> = [1, 2, 4, 8]
        .into_iter()
        .filter(|base| mask & base != 0)
        .collect();
    if bases.is_empty() {
        return letter(rng);
    }
    bases[rng.below(bases.len())]
}
