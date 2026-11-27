mod common;

use _kernels::indel::{self, BLOCK};
use _kernels::letters::Letters;
use common::{SplitMix64, compat, copy_with_edits, letter, plant, sequence};

const BLOCKS: [usize; 7] = [1, 2, 3, 5, 8, 64, 301];

struct Case {
    id: usize,
    seq: Vec<u8>,
    letters: Letters,
    max_mm: usize,
    block: usize,
}

/// 20,000 cases from a fixed seed: m from 2 to 64, max_mm from 1 to m - 1, a sequence of 0 to 300
/// letters with 0 to 2 planted copies, a random pair of fixed flags, and a block from BLOCKS.
fn cases() -> impl Iterator<Item = Case> {
    let mut rng = SplitMix64::new(4);
    (0..20_000).map(move |id| {
        let m = 2 + rng.below(63);
        let pattern: Vec<u8> = (0..m).map(|_| letter(&mut rng)).collect();
        let len = rng.below(301);
        let mut seq = sequence(&mut rng, len);
        for _ in 0..rng.below(3) {
            let edits = rng.below(4);
            plant(&mut rng, &mut seq, &pattern, edits);
        }
        let compat = compat(rng.one_in(2), rng.one_in(2));
        Case {
            id,
            seq,
            letters: Letters::new(&pattern, &compat),
            max_mm: 1 + rng.below(m - 1),
            block: BLOCKS[rng.below(BLOCKS.len())],
        }
    })
}

fn bits_and_scalar(
    seq: &[u8],
    letters: &Letters,
    max_mm: usize,
    block: usize,
) -> [Vec<(i64, i64)>; 2] {
    let (mut bits, mut scalar) = (Vec::new(), Vec::new());
    indel::bits(seq, letters, max_mm, block, &mut |s, e| bits.push((s, e)));
    indel::scalar(seq, letters, max_mm, &mut |s, e| scalar.push((s, e)));
    [bits, scalar]
}

#[test]
fn bits_equal_scalar_on_random_input() {
    for case in cases() {
        let [bits, scalar] = bits_and_scalar(&case.seq, &case.letters, case.max_mm, case.block);
        assert_eq!(bits, scalar, "case {}", case.id);
    }
}

#[test]
fn filter_marks_each_accepted_position() {
    for case in cases() {
        let candidates = indel::candidates(&case.seq, &case.letters, case.max_mm, case.block);
        for j in indel::accepted(&case.seq, &case.letters, case.max_mm) {
            assert!(
                candidates.binary_search(&j).is_ok(),
                "case {} position {j}",
                case.id
            );
        }
    }
}

#[test]
fn bits_equal_scalar_on_a_sequence_longer_than_one_block() {
    let mut rng = SplitMix64::new(5);
    let base = |rng: &mut SplitMix64| [1, 2, 4, 8][rng.below(4)];
    let pattern: Vec<u8> = (0..21).map(|_| base(&mut rng)).collect();
    let mut seq: Vec<u8> = (0..3 * BLOCK + 17).map(|_| base(&mut rng)).collect();
    // The copy at BLOCK writes over the copy at BLOCK - 1, and the copy at 3 * BLOCK keeps only
    // 17 of its 21 letters. Thus that copy has no edit, and it matches with 4 deletions.
    let copies = [
        (BLOCK - 21, 1),
        (BLOCK - 1, 4),
        (BLOCK, 2),
        (2 * BLOCK - 10, 3),
        (3 * BLOCK, 0),
    ];
    for (start, edits) in copies {
        let copy = copy_with_edits(&mut rng, &pattern, edits);
        for (slot, letter) in seq[start..].iter_mut().zip(copy) {
            *slot = letter;
        }
    }
    let letters = Letters::new(&pattern, &compat(true, true));
    for max_mm in 1..=4 {
        let [bits, scalar] = bits_and_scalar(&seq, &letters, max_mm, BLOCK);
        assert_eq!(bits, scalar, "max_mm {max_mm}");
    }
    let [_, scalar] = bits_and_scalar(&seq, &letters, 4, BLOCK);
    let starts_at = |j: usize| scalar.iter().any(|&(start, _)| start == one_based(j));
    let crosses = |j: usize| {
        let limit = one_based(j);
        scalar
            .iter()
            .any(|&(start, end)| start < limit && end >= limit)
    };
    assert!(starts_at(BLOCK), "{scalar:?}");
    assert!(crosses(2 * BLOCK), "{scalar:?}");
    assert!(starts_at(3 * BLOCK), "{scalar:?}");
}

fn one_based(j: usize) -> i64 {
    i64::try_from(j + 1).expect("a position of the test sequence fits in i64")
}
