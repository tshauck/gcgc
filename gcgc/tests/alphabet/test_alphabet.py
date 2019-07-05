# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

import unittest

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import pytest

from gcgc import alphabet
from gcgc.alphabet.utils import biopython_alphabet_to_gcgc_alphabet
from gcgc.encoded_seq import EncodedSeq
from gcgc.exceptions import GCGCAlphabetLetterEncodingException


@pytest.mark.parametrize(
    "biopython_class,gcgc_class",
    [
        (IUPAC.IUPACUnambiguousDNA, alphabet.IUPACUnambiguousDNAEncoding),
        (IUPAC.IUPACAmbiguousDNA, alphabet.IUPACAmbiguousDNAEncoding),
        (IUPAC.IUPACUnambiguousRNA, alphabet.IUPACUnambiguousRNAEncoding),
        (IUPAC.IUPACAmbiguousRNA, alphabet.IUPACAmbiguousRNAEncoding),
        (IUPAC.ExtendedIUPACDNA, alphabet.ExtendedIUPACDNAEncoding),
        (IUPAC.ExtendedIUPACProtein, alphabet.ExtendedIUPACProteinEncoding),
        (IUPAC.IUPACProtein, alphabet.IUPACProteinEncoding),
    ],
)
def test_biopython_alphabet_to_gcgc_alphabet(biopython_class, gcgc_class):
    biopython_instance = biopython_class()

    klass = biopython_alphabet_to_gcgc_alphabet(biopython_instance)
    assert isinstance(klass, biopython_class)


@pytest.mark.parametrize(
    "biopython_class,gcgc_class",
    [
        (IUPAC.IUPACUnambiguousDNA, alphabet.IUPACUnambiguousDNAEncoding),
        (IUPAC.IUPACAmbiguousDNA, alphabet.IUPACAmbiguousDNAEncoding),
        (IUPAC.IUPACUnambiguousRNA, alphabet.IUPACUnambiguousRNAEncoding),
        (IUPAC.IUPACAmbiguousRNA, alphabet.IUPACAmbiguousRNAEncoding),
        (IUPAC.ExtendedIUPACDNA, alphabet.ExtendedIUPACDNAEncoding),
        (IUPAC.ExtendedIUPACProtein, alphabet.ExtendedIUPACProteinEncoding),
        (IUPAC.IUPACProtein, alphabet.IUPACProteinEncoding),
    ],
)
def test_biopython_from_seq(biopython_class, gcgc_class):
    biopython_instance = biopython_class()

    es = EncodedSeq.from_seq(Seq("ATCG", biopython_instance))
    assert isinstance(es.alphabet, gcgc_class)


@pytest.mark.parametrize("kmer_size,start,expected_len", [(1, True, 7), (2, False, 18)])
def test_len(kmer_size, start, expected_len):
    dna = alphabet.IUPACUnambiguousDNAEncoding(kmer_size=kmer_size, start_token=start)
    assert len(dna) == expected_len


class TestAlphabet(unittest.TestCase):
    def test_decoding_index(self):
        dna = alphabet.IUPACUnambiguousDNAEncoding()
        self.assertEqual(dna.decode_token(0), dna.decoding_index[0])

    def test_encoding_index(self):
        dna = alphabet.IUPACUnambiguousDNAEncoding()

        self.assertEqual(dna.encode_token("A"), dna.encoding_index["A"])

    def test_raise_integer_encoding_exception(self):

        dna = alphabet.IUPACUnambiguousDNAEncoding()

        with self.assertRaises(GCGCAlphabetLetterEncodingException):
            dna.integer_encode("Z")

    def test_integer_decode(self):

        dna = alphabet.IUPACUnambiguousDNAEncoding()

        code = [0, 1, 2, 3]
        expected = "".join(dna.decode_token(s) for s in code)
        actual = dna.integer_decode(code)

        self.assertEqual(expected, actual)

    def test_kmer_tokens_size(self):
        dna = alphabet.IUPACUnambiguousDNAEncoding(kmer_size=2)
        n_kmers = len(dna.kmers)
        n_kmers_and_tokens = len(dna.kmers_and_tokens)

        self.assertEqual(n_kmers, 4 ** 2)
        self.assertEqual(n_kmers_and_tokens, 4 ** 2 + 3)


@pytest.mark.parametrize(
    "seq,kmer_size,expected_kmer",
    [
        ("ATCG", 2, ["AT", "TC", "CG"]),
        ("ATCGAT", 3, ["ATC", "TCG", "CGA", "GAT"]),
        ("ATCG", 1, ["A", "T", "C", "G"]),
    ],
)
def test_kmer_encoding(seq, kmer_size, expected_kmer):
    """Test the kemrs are encoded as expected."""
    dna = alphabet.IUPACUnambiguousDNAEncoding(kmer_size=kmer_size)
    expected_integers = [dna.encode_token(t) for t in expected_kmer]

    actual = dna.integer_encode(seq)

    assert expected_integers == actual


def test_special_token_integer_encoding():
    """Test the special characters encode to integers correctly."""
    dna = alphabet.IUPACUnambiguousDNAEncoding()

    assert dna.encoded_start == dna.encode_token(dna.START)
    assert dna.encoded_end == dna.encode_token(dna.END)

    assert dna.encoded_padding == dna.encode_token(dna.PADDING)
    assert dna.encoded_padding == 0


@pytest.mark.parametrize(
    "start_token,end_token,expected_tokens",
    [(False, True, "|<"), (False, False, "|"), (True, False, "|>")],
)
def test_alphabet_configuration(start_token, end_token, expected_tokens):
    """Test that we can selectively use start and end tokens."""
    dna = alphabet.IUPACUnambiguousDNAEncoding(start_token=start_token, end_token=end_token)
    assert dna.tokens == expected_tokens
