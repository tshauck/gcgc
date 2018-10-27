# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest

from Bio.Alphabet import IUPAC

from gcgc.alphabet.utils import biopython_alphabet_to_gcgc_alphabet
from gcgc.alphabet import IUPACDnaAlphabetEncoding


class TestUtils(unittest.TestCase):
    def test_biopython_alphabet_to_gcgc_alphabet(self):
        dna = IUPAC.IUPACUnambiguousDNA()

        klass = biopython_alphabet_to_gcgc_alphabet(dna)
        self.assertIsInstance(klass, IUPACDnaAlphabetEncoding)
