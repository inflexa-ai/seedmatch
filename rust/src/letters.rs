//! The letter table: which subject masks each pattern letter matches.

pub const MASKS: usize = 16;
pub const WORD: usize = 64;

/// `compat[p][s]` is true if the pattern mask `p` matches the subject mask `s`.
pub type Compat = [[bool; MASKS]; MASKS];

pub struct Letters {
    ok: Vec<[bool; MASKS]>,
    first: [usize; MASKS],
}

impl Letters {
    pub fn new(pattern: &[u8], compat: &Compat) -> Self {
        let ok: Vec<[bool; MASKS]> = pattern
            .iter()
            .map(|&mask| compat[usize::from(mask & 15)])
            .collect();
        let mut first = [ok.len(); MASKS];
        for (s, slot) in first.iter_mut().enumerate() {
            if let Some(a) = ok.iter().position(|row| row[s]) {
                *slot = a;
            }
        }
        Self { ok, first }
    }

    pub fn len(&self) -> usize {
        self.ok.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ok.is_empty()
    }

    /// Return true if pattern letter `a` matches the subject code.
    pub fn ok(&self, a: usize, code: u8) -> bool {
        self.ok[a][usize::from(code & 15)]
    }

    /// Return the first pattern letter that matches the subject code, or `len()` if none does.
    pub fn first(&self, code: u8) -> usize {
        self.first[usize::from(code & 15)]
    }

    /// Return one word for each subject mask s: bit a is 1 if pattern letter a matches s. The
    /// letters after the first `WORD` letters are not in the words.
    pub fn words(&self) -> [u64; MASKS] {
        let mut words = [0; MASKS];
        for (a, row) in self.ok.iter().take(WORD).enumerate() {
            for (word, &ok) in words.iter_mut().zip(row) {
                if ok {
                    *word |= 1 << a;
                }
            }
        }
        words
    }
}
