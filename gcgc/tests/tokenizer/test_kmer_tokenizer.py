# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

import pytest

from gcgc import KmerTokenizer


@pytest.mark.parametrize(
    "seq,tokenizer,expected_tokens,expected_encoding",
    [
        (
            "ATCG",
            KmerTokenizer.bare_tokenizer(alphabet="ATCG", kmer_length=1, kmer_stride=1),
            list("ATCG"),
            [0, 1, 2, 3],
        ),
        (
            "ATCG",
            KmerTokenizer.bare_tokenizer(
                alphabet="ATCG", max_length=2, kmer_length=1, kmer_stride=1
            ),
            list("AT"),
            [0, 1],
        ),
        (
            "ATCG",
            KmerTokenizer.bare_tokenizer(
                alphabet="ATCG", min_length=4, pad_token="|", kmer_length=2, kmer_stride=1
            ),
            ["AT", "TC", "CG", "|"],
            [2, 7, 12, 0],
        ),
        (
            "AATT",
            KmerTokenizer.bare_tokenizer(
                alphabet="ATCG", pad_token="|", conform_length=5, kmer_length=2, kmer_stride=2
            ),
            ["AA", "TT", "|", "|", "|"],
            [1, 6, 0, 0, 0],
        ),
        (
            "AATT",
            KmerTokenizer.bare_tokenizer(
                alphabet="ATCG",
                pad_token="|",
                conform_length=5,
                kmer_length=2,
                kmer_stride=2,
                pad_at_end=False,
            ),
            ["|", "|", "|", "AA", "TT"],
            [0, 0, 0, 1, 6],
        ),
        (
            "ATCGATCGATCGATCG",
            KmerTokenizer.bare_tokenizer(
                alphabet="ATCG", max_length=2, kmer_length=2, kmer_stride=1
            ),
            ["AT", "TC"],
            [1, 6],
        ),
        (
            "ATCGATCG",
            KmerTokenizer(
                alphabet="ATCG",
                pad_token="|",
                bos_token=">",
                eos_token="<",
                unk_token=None,
                mask_token=None,
                max_length=10,
                kmer_length=2,
                kmer_stride=1,
            ),
            [">", "AT", "TC", "CG", "GA", "AT", "TC", "CG", "<"],
            [1, 4, 9, 14, 15, 4, 9, 14, 2],
        ),
    ],
)
def test_kmer_tokenization(seq, tokenizer, expected_tokens, expected_encoding):
    """Test the kmers are encoded as expected."""
    actual_tokens = tokenizer.tokenize(seq)
    assert actual_tokens == expected_tokens

    actual_encoded = tokenizer.encode(seq)
    assert actual_encoded == expected_encoding


@pytest.mark.parametrize(
    "token_ids, tokenizer, expected_mask",
    [
        (
            [0, 1, 2, 3, 4],
            KmerTokenizer(
                alphabet="ATCG",
                kmer_length=1,
                kmer_stride=1,
                pad_token="|",
                bos_token=None,
                eos_token=None,
                unk_token=None,
                mask_token=None,
            ),
            [1, 0, 0, 0, 0],
        ),
    ],
)
def test_get_special_tokens_mask(token_ids, tokenizer, expected_mask):
    """Test the special tokens mask returns the proper mask for a given input list of token ids."""
    actual_mask = tokenizer.get_special_tokens_mask(token_ids)
    assert actual_mask == expected_mask
