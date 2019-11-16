# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

import pytest

from gcgc.tokenizer import SequenceTokenizer, SequenceTokenizerSpec


@pytest.mark.parametrize(
    "seq,spec,expected_encoding",
    [
        ("ATCG", SequenceTokenizerSpec(max_length=10, kmer_size=1, kmer_step_size=1), list("ATCG")),
        (
            "ATCG",
            SequenceTokenizerSpec(max_length=10, kmer_size=2, kmer_step_size=1),
            ["AT", "TC", "CG"],
        ),
        ("ATCG", SequenceTokenizerSpec(max_length=10, kmer_size=2, kmer_step_size=2), ["AT", "CG"]),
    ],
)
def test_kmer_tokenization(seq, spec, expected_encoding):
    """Test the kemrs are encoded as expected."""
    tokenizer = SequenceTokenizer(spec)
    actual = tokenizer(seq)

    assert actual == expected_encoding
