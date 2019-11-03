# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

import pytest

from gcgc.tokenizer import SequenceTokenizer


@pytest.mark.parametrize(
    "seq,kmer_size,kmer_step_size,expected_encoding",
    [
        ("ATCG", 1, 1, list("ATCG")),
        ("ATCG", 2, 1, ["AT", "TC", "CG"]),
        ("ATCG", 2, 2, ["AT", "CG"]),
    ],
)
def test_kmer_tokenization(seq, kmer_size, kmer_step_size, expected_encoding):
    """Test the kemrs are encoded as expected."""
    tokenizer = SequenceTokenizer(kmer_size, kmer_step_size)
    actual = tokenizer(seq)

    assert actual == expected_encoding
