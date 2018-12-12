# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest

from Bio.Alphabet import IUPAC

from gcgc.alphabet import IUPACUnambiguousDNAEncoding
from gcgc.alphabet.utils import biopython_alphabet_to_gcgc_alphabet
from gcgc.exceptions import GCGCAlphabetLetterEncodingException


class TestAlphabet(unittest.TestCase):
    def test_biopython_alphabet_to_gcgc_alphabet(self):
        dna = IUPAC.IUPACUnambiguousDNA()

        klass = biopython_alphabet_to_gcgc_alphabet(dna)
        self.assertIsInstance(klass, IUPACUnambiguousDNAEncoding)

    def test_len(self):
        dna = IUPACUnambiguousDNAEncoding()
        self.assertEqual(len(dna), len(dna.letters_and_tokens))

    def test_decoding_index(self):
        dna = IUPACUnambiguousDNAEncoding()

        self.assertEqual(dna.decode_token(0), dna.decoding_index[0])

    def test_encoding_index(self):
        dna = IUPACUnambiguousDNAEncoding()

        self.assertEqual(dna.encode_token("A"), dna.encoding_index["A"])

    def test_raise_integer_encoding_exception(self):

        dna = IUPACUnambiguousDNAEncoding()

        with self.assertRaises(GCGCAlphabetLetterEncodingException):
            dna.integer_encode("Z")

    def test_integer_decode(self):

        dna = IUPACUnambiguousDNAEncoding()

        code = [0, 1, 2, 3]
        expected = "".join(dna.decode_token(s) for s in code)
        actual = dna.integer_decode(code)

        self.assertEqual(expected, actual)
