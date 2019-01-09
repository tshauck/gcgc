# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest

from Bio import AlignIO
from Bio.Alphabet import IUPAC

from gcgc.encoded_seq import EncodedSeq
from gcgc.tests.fixtures import PF01152_PATH


class TestAlignment(unittest.TestCase):
    def test_read_alignment_file(self):

        protein_alphabet = IUPAC.ExtendedIUPACProtein()

        gaps = {".", "-"}

        # Times two because add_lower_case_for_inserts=True, 3 because start, end, padding.
        expected_dim = 2 * (len(protein_alphabet.letters)) + len(gaps) + 3

        with open(PF01152_PATH) as f:

            alignment = AlignIO.read(f, "stockholm", alphabet=protein_alphabet)
            for f in alignment:

                seq = EncodedSeq.from_seq(
                    f.seq, gap_characters=gaps, add_lower_case_for_inserts=True
                )

                one_hot_dim = len(seq.one_hot_encoded[0])
                self.assertEqual(one_hot_dim, expected_dim)
