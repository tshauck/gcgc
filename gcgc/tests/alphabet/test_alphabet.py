# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest

import pytest
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

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


class TestAlphabet(unittest.TestCase):
    def test_len(self):
        dna = alphabet.IUPACUnambiguousDNAEncoding()
        self.assertEqual(len(dna), len(dna.letters_and_tokens))

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
