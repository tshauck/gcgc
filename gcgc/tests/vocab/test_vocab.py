# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Test the vocabulary object."""

from typing import NamedTuple
from typing import Set

from gcgc import KmerTokenizer
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


def test_vocab_values():
    """Test that for the alphabet the correct vocab is generated."""

    class Case(NamedTuple):
        kmer_tokenzier: KmerTokenizer
        expected: Set

    test_cases = [
        Case(KmerTokenizer.bare_tokenizer(alphabet="ATGC"), {"A", "T", "C", "G"}),
        Case(KmerTokenizer.bare_tokenizer(alphabet="ATGCN"), {"A", "T", "C", "G", "N"}),
        Case(KmerTokenizer.bare_tokenizer(alphabet="extended_dna"), set("GATCBDSW")),
        Case(KmerTokenizer.bare_tokenizer(alphabet="ambiguous_dna"), set("GATCRYWSMKHBVDN")),
    ]

    for case in test_cases:
        assert set(case.kmer_tokenzier.vocab.stoi.keys()) == case.expected, case.kmer_tokenzier
