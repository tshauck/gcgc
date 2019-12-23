# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

import pytest

from gcgc.tokenizer import KmerTokenizer, KmerTokenizerSettings


@pytest.mark.parametrize(
    "seq,settings,expected_tokens,expected_encoding",
    [
        (
            "ATCG",
            KmerTokenizerSettings(alphabet="ATCG", kmer_length=1, kmer_stride=1),
            list("ATCG"),
            [0, 1, 2, 3],
        ),
        (
            "ATCG",
            KmerTokenizerSettings(alphabet="ATCG", max_length=2, kmer_length=1, kmer_stride=1),
            list("AT"),
            [0, 1],
        ),
        (
            "ATCG",
            KmerTokenizerSettings(alphabet="ATCG", kmer_length=2, kmer_stride=1),
            ["AT", "TC", "CG"],
            [1, 6, 11],
        ),
        (
            "AATT",
            KmerTokenizerSettings(alphabet="ATCG", max_length=10, kmer_length=2, kmer_stride=2),
            ["AA", "TT"],
            [0, 5],
        ),
        (
            "ATCGATCGATCGATCG",
            KmerTokenizerSettings(alphabet="ATCG", max_length=2, kmer_length=2, kmer_stride=1),
            ["AT", "TC"],
            [1, 6],
        ),
        (
            "ATCGATCG",
            KmerTokenizerSettings(
                alphabet="ATCG",
                bos_token=">",
                eos_token="<",
                max_length=10,
                kmer_length=2,
                kmer_stride=1,
            ),
            [">", "AT", "TC", "CG", "GA", "AT", "TC", "CG", "<"],
            [0, 3, 8, 13, 14, 3, 8, 13, 1],
        ),
    ],
)
def test_kmer_tokenization(seq, settings, expected_tokens, expected_encoding):
    """Test the kemrs are encoded as expected."""
    tokenizer = KmerTokenizer(settings)
    actual_tokens = tokenizer.encode_as_tokens(seq)
    assert actual_tokens == expected_tokens

    actual_encoded = tokenizer.encode_as_ids(seq)
    assert actual_encoded == expected_encoding
