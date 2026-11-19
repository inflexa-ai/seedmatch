from hypothesis import strategies as st


def dna(letters, *, max_size, min_size=0):
    # st.text keeps each letter of its alphabet once, so a repeated letter would not add weight.
    # st.lists gives mostly short lists, so the length comes first, from a uniform range.
    return st.integers(min_size, max_size).flatmap(
        lambda size: st.lists(st.sampled_from(letters), min_size=size, max_size=size).map("".join)
    )


@st.composite
def subjects_for(draw, pattern, max_size):
    """A random subject. In half of the cases, it holds a copy of the pattern."""
    subject = draw(dna("ACGTN" * 6 + "WSK", max_size=max_size))
    if draw(st.booleans()):
        at = draw(st.integers(0, len(subject)))
        subject = subject[:at] + pattern + subject[at:]
    return subject
