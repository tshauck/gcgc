# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

import pytest

from gcgc.tokenizer import KmerTokenizer
from gcgc.tokenizer import KmerTokenizerSettings


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
            KmerTokenizerSettings(
                alphabet="ATCG", min_length=4, pad_token="|", kmer_length=2, kmer_stride=1
            ),
            ["AT", "TC", "CG", "|"],
            [2, 7, 12, 0],
        ),
        (
            "AATT",
            KmerTokenizerSettings(
                alphabet="ATCG", pad_token="|", conform_length=5, kmer_length=2, kmer_stride=2
            ),
            ["AA", "TT", "|", "|", "|"],
            [1, 6, 0, 0, 0],
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
            [1, 3, 8, 13, 14, 3, 8, 13, 2],
        ),
    ],
)
def test_kmer_tokenization(seq, settings, expected_tokens, expected_encoding):
    """Test the kmers are encoded as expected."""
    tokenizer = KmerTokenizer(settings)
    actual_tokens = tokenizer.encode_as_tokens(seq)
    assert actual_tokens == expected_tokens

    actual_encoded = tokenizer.encode(seq)
    assert actual_encoded == expected_encoding


@pytest.mark.parametrize(
    "token_ids, settings, expected_mask",
    [
        (
            [0, 1, 2, 3],
            KmerTokenizerSettings(alphabet="ATCG", kmer_length=1, kmer_stride=1),
            [0, 0, 0, 0],
        ),
        (
            [0, 1, 2, 3, 4],
            KmerTokenizerSettings(alphabet="ATCG", pad_token="|", kmer_length=1, kmer_stride=1),
            [1, 0, 0, 0, 0],
        ),
        (
            [0, 1, 2, 3, 4],
            KmerTokenizerSettings(
                alphabet="ATCG", bos_token_id=0, bos_token=">", kmer_length=1, kmer_stride=1
            ),
            [1, 0, 0, 0, 0],
        ),
    ],
)
def test_get_special_tokens_mask(token_ids, settings, expected_mask):
    """Test the special tokens mask returns the proper mask for a given input list of token ids."""

    tokenizer = KmerTokenizer(settings)
    actual_mask = tokenizer.get_special_tokens_mask(token_ids)
    assert actual_mask == expected_mask
