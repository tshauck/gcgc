# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Test the vocabulary object."""

from gcgc.vocab import Vocab


def test_vocab():
    """Test the vocabulary object works as expected."""
    vocab = Vocab()

    vocab.add_item("A")

    assert "A" in vocab
    assert vocab["A"] == 0

    vocab.add_items(["T", "C"])

    assert "T" in vocab
    assert "C" in vocab

    assert vocab.stoi == {"A": 0, "T": 1, "C": 2}
