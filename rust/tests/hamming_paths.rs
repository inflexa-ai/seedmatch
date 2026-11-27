mod common;

use _kernels::hamming;
use _kernels::letters::Letters;
use common::{SplitMix64, compat, letter, plant, sequence};

struct Case {
    id: usize,
    seq: Vec<u8>,
    letters: Letters,
    max_mm: u64,
    min_mm: u64,
}

/// 20,000 cases from a fixed seed. Each case: m from 1 to 64, a sequence of 0 to 300 letters with
/// 0 to 2 planted copies, a random pair of fixed flags, max_mm from 0 to m + 2, and min_mm from 0
/// to max_mm. One case in 50 uses max_mm = min_mm = u64::MAX.
fn cases() -> impl Iterator<Item = Case> {
    let mut rng = SplitMix64::new(3);
    (0..20_000).map(move |id| {
        let m = 1 + rng.below(64);
        let pattern: Vec<u8> = (0..m).map(|_| letter(&mut rng)).collect();
        let len = rng.below(301);
        let mut seq = sequence(&mut rng, len);
        for _ in 0..rng.below(3) {
            let edits = rng.below(4);
            plant(&mut rng, &mut seq, &pattern, edits);
        }
        let compat = compat(rng.one_in(2), rng.one_in(2));
        let (max_mm, min_mm) = if rng.one_in(50) {
            (u64::MAX, u64::MAX)
        } else {
            let max_mm = rng.below(m + 3);
            (max_mm as u64, rng.below(max_mm + 1) as u64)
        };
        Case {
            id,
            seq,
            letters: Letters::new(&pattern, &compat),
            max_mm,
            min_mm,
        }
    })
}

#[test]
fn bits_equal_scalar_on_random_input() {
    for case in cases() {
        let (mut bits, mut scalar) = (Vec::new(), Vec::new());
        hamming::bits(
            &case.seq,
            &case.letters,
            case.max_mm,
            case.min_mm,
            &mut |s| bits.push(s),
        );
        hamming::scalar(
            &case.seq,
            &case.letters,
            case.max_mm,
            case.min_mm,
            &mut |s| scalar.push(s),
        );
        assert_eq!(bits, scalar, "case {}", case.id);
    }
}
