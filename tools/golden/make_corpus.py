"""Write the search cases that tools/golden/R/search_corpus.R gives to Biostrings.

`corpus_cases.tsv` holds short random cases with many edge cases. `seed_corpus_cases.tsv` holds
patterns with the lengths of the seeds and longer subjects with edited copies of the pattern.
"""

import random
from pathlib import Path

SEED = 20261119
CASES = 400
MAX_PATTERN = 12
MAX_SUBJECT = 60
MAX_EDITS = 3
BASES = "ACGT"
AMBIGUOUS = "MRWSYKVHDB"
N_RUN_RATE = 0.03
END_EDIT_RATE = 0.3
PATTERN_FROM_REPEAT_RATE = 0.6
PART_OF_COPY_RATE = 0.5
SEED_CORPUS_SEED = 20261120
SEED_CORPUS_CASES = 1000
SEED_LENGTHS = ((6, 8), (15, 23))
SEED_SUBJECT = (20, 150)
SEED_EDITS = (1, 4)
MAX_COPIES = 3
END_COPY_RATE = 0.2
GOLDEN = Path(__file__).resolve().parents[2] / "tests" / "golden" / "search"


def random_letters(rng: random.Random, length: int, *, n_run: int, ambiguous: float) -> str:
    """Mostly ACGT, with runs of up to `n_run` N letters and some other ambiguity codes."""
    letters: list[str] = []
    while len(letters) < length:
        roll = rng.random()
        if roll < N_RUN_RATE:
            letters.extend("N" * rng.randint(1, n_run))
        elif roll < N_RUN_RATE + ambiguous:
            letters.append(rng.choice(AMBIGUOUS))
        else:
            letters.append(rng.choice(BASES))
    return "".join(letters[:length])


def pattern_letters(rng: random.Random, length: int) -> str:
    return random_letters(rng, length, n_run=3, ambiguous=0.05)


def subject_letters(rng: random.Random, length: int) -> str:
    return random_letters(rng, length, n_run=8, ambiguous=0.03)


def edited(rng: random.Random, seq: str, edits: int) -> str:
    """Apply substitutions, insertions, and deletions. Some edits go at the first or last place."""
    letters = list(seq)
    for _ in range(edits):
        kind = rng.choice(("substitute", "insert", "delete"))
        if kind == "delete" and len(letters) <= 1:
            kind = "substitute"
        size = len(letters) + 1 if kind == "insert" else len(letters)
        place = rng.choice((0, size - 1)) if rng.random() < END_EDIT_RATE else rng.randrange(size)
        if kind == "substitute":
            letters[place] = subject_letters(rng, 1)
        elif kind == "insert":
            letters.insert(place, subject_letters(rng, 1))
        else:
            del letters[place]
    return "".join(letters)


def copy_of(rng: random.Random, pattern: str) -> str:
    return edited(rng, pattern, rng.randint(0, MAX_EDITS))


def split_flanks(rng: random.Random, budget: int) -> tuple[int, int]:
    total = rng.randint(0, budget)
    left = rng.randint(0, total)
    return left, total - left


def planted_once(rng: random.Random, pattern: str) -> str:
    """A copy inside the subject, at one end, or cut at one end."""
    copy = copy_of(rng, pattern)
    place = rng.choice(("inside", "first", "last", "cut_first", "cut_last"))
    cut = rng.randint(1, max(1, len(copy) // 3)) if len(copy) > 1 else 0
    if place == "cut_first":
        copy = copy[cut:]
    elif place == "cut_last":
        copy = copy[: len(copy) - cut]
    left, right = split_flanks(rng, MAX_SUBJECT - len(copy))
    if place in ("first", "cut_first"):
        left = 0
    if place in ("last", "cut_last"):
        right = 0
    return subject_letters(rng, left) + copy + subject_letters(rng, right)


def planted_twice(rng: random.Random, pattern: str) -> str:
    """Two copies with a short gap, so that near matches compete."""
    first, second = copy_of(rng, pattern), copy_of(rng, pattern)
    gap = subject_letters(rng, rng.randint(0, 4))
    middle = first + gap + second
    left, right = split_flanks(rng, MAX_SUBJECT - len(middle))
    return subject_letters(rng, left) + middle + subject_letters(rng, right)


def repeated(rng: random.Random, length: int) -> tuple[str, str]:
    """A subject of a repeated short unit. Most patterns come from the same repeat."""
    unit = random_letters(rng, rng.randint(1, 4), n_run=1, ambiguous=0.05)
    subject = (unit * MAX_SUBJECT)[: rng.randint(0, MAX_SUBJECT)]
    if rng.random() < PATTERN_FROM_REPEAT_RATE:
        offset = rng.randrange(len(unit))
        part = (unit * MAX_PATTERN)[offset : offset + length]
        return edited(rng, part, rng.randint(0, 1))[:MAX_PATTERN], subject
    return pattern_letters(rng, length), subject


def shorter_than_pattern(rng: random.Random, pattern: str) -> str:
    """An empty subject, or a part of a copy of the pattern, or random letters."""
    length = rng.randint(0, len(pattern) - 1)
    if length > 0 and rng.random() < PART_OF_COPY_RATE:
        copy = copy_of(rng, pattern)
        start = rng.randint(0, max(0, len(copy) - length))
        return copy[start : start + length]
    return subject_letters(rng, length)


def make_case(rng: random.Random) -> tuple[str, str]:
    length = rng.randint(1, MAX_PATTERN)
    kind = rng.choices(("once", "twice", "repeat", "random", "short"), (55, 12, 12, 11, 10))[0]
    if kind == "repeat":
        return repeated(rng, length)
    pattern = pattern_letters(rng, length)
    if kind == "once":
        return pattern, planted_once(rng, pattern)
    if kind == "twice":
        return pattern, planted_twice(rng, pattern)
    if kind == "short":
        return pattern, shorter_than_pattern(rng, pattern)
    return pattern, subject_letters(rng, rng.randint(0, MAX_SUBJECT))


def seed_case(rng: random.Random) -> tuple[str, str]:
    """A pattern with the length of a seed, and a subject with edited copies of it."""
    low, high = rng.choice(SEED_LENGTHS)
    pattern = pattern_letters(rng, rng.randint(low, high))
    copies = [
        edited(rng, pattern, rng.randint(*SEED_EDITS)) for _ in range(rng.randint(1, MAX_COPIES))
    ]
    length = rng.randint(*SEED_SUBJECT)
    background = subject_letters(rng, max(0, length - sum(len(copy) for copy in copies)))
    cuts = sorted(rng.randint(0, len(background)) for _ in copies)
    if rng.random() < END_COPY_RATE:
        cuts[0] = 0
    if rng.random() < END_COPY_RATE:
        cuts[-1] = len(background)
    parts, last = [], 0
    for cut, copy in zip(cuts, copies, strict=True):
        parts += [background[last:cut], copy]
        last = cut
    parts.append(background[last:])
    return pattern, "".join(parts)


def seeded(seed: int) -> random.Random:
    # The seed must give the same corpus again, so the generator is seeded and not secure.
    return random.Random(seed)  # noqa: S311


def check_bounds(cases: list[tuple[str, str]], patterns: range, subjects: range) -> None:
    for case_id, (pattern, subject) in enumerate(cases, start=1):
        if len(pattern) not in patterns or len(subject) not in subjects:
            message = f"case {case_id}: pattern {pattern!r} or subject {subject!r} out of bounds"
            raise ValueError(message)


def random_cases() -> list[tuple[str, str]]:
    rng = seeded(SEED)
    cases = [make_case(rng) for _ in range(CASES)]
    check_bounds(cases, range(1, MAX_PATTERN + 1), range(MAX_SUBJECT + 1))
    return cases


def seed_cases() -> list[tuple[str, str]]:
    rng = seeded(SEED_CORPUS_SEED)
    cases = [seed_case(rng) for _ in range(SEED_CORPUS_CASES)]
    patterns = range(SEED_LENGTHS[0][0], SEED_LENGTHS[-1][1] + 1)
    check_bounds(cases, patterns, range(SEED_SUBJECT[0], SEED_SUBJECT[1] + 1))
    return cases


def write_cases(path: Path, seed: int, cases: list[tuple[str, str]]) -> None:
    lines = [
        f"# generated-by: tools/golden/make_corpus.py with random.Random({seed})",
        "case_id\tpattern\tsubject",
    ]
    lines += [
        f"{case_id}\t{pattern}\t{subject}" for case_id, (pattern, subject) in enumerate(cases, 1)
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    for name, seed, cases in (
        ("corpus", SEED, random_cases()),
        ("seed_corpus", SEED_CORPUS_SEED, seed_cases()),
    ):
        path = GOLDEN / f"{name}_cases.tsv"
        write_cases(path, seed, cases)
        print(f"wrote {path}")
