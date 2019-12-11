# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

import pytest

from gcgc.tokenizer import SequenceTokenizer, SequenceTokenizerSpec


@pytest.mark.parametrize(
    "seq,spec,expected_tokens,expected_encoding",
    [
        (
            "ATCG",
            SequenceTokenizerSpec(alphabet="ATCG", max_length=10, kmer_size=1, kmer_step_size=1),
            list("ATCG"),
            [0, 1, 2, 3],
        ),
        (
            "ATCG",
            SequenceTokenizerSpec(alphabet="ATCG", kmer_size=1, kmer_step_size=1),
            list("ATCG"),
            [0, 1, 2, 3],
        ),
        (
            "ATCG",
            SequenceTokenizerSpec(alphabet="ATCG", max_length=2, kmer_size=1, kmer_step_size=1),
            list("AT"),
            [0, 1],
        ),
        (
            "ATCG",
            SequenceTokenizerSpec(alphabet="ATCG", max_length=10, kmer_size=2, kmer_step_size=1),
            ["AT", "TC", "CG"],
            [1, 6, 11],
        ),
        (
            "AATT",
            SequenceTokenizerSpec(alphabet="ATCG", max_length=10, kmer_size=2, kmer_step_size=2),
            ["AA", "TT"],
            [0, 5],
        ),
        (
            "ATCGATCGATCGATCG",
            SequenceTokenizerSpec(alphabet="ATCG", max_length=2, kmer_size=2, kmer_step_size=1),
            ["AT", "TC"],
            [1, 6],
        ),
        (
            "ATCGATCG",
            SequenceTokenizerSpec(
                alphabet="ATCG",
                bos_token=">",
                eos_token="<",
                max_length=10,
                kmer_size=2,
                kmer_step_size=1,
            ),
            [">", "AT", "TC", "CG", "GA", "AT", "TC", "CG", "<"],
            [0, 3, 8, 13, 14, 3, 8, 13, 1],
        ),
    ],
)
def test_kmer_tokenization(seq, spec, expected_tokens, expected_encoding):
    """Test the kemrs are encoded as expected."""
    tokenizer = SequenceTokenizer(spec)
    actual_tokens = tokenizer(seq)
    assert actual_tokens == expected_tokens

    actual_encoded = tokenizer.encode(seq)
    assert actual_encoded == expected_encoding
