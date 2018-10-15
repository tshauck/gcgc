"""Unit test package for gcgc."""

import unittest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np

from gcgc.alphabet.iupac import UnambiguousDnaAlphabet
from gcgc.alphabet.iupac import ExtendedProteinAlphabet
from gcgc.seq_record import EncodedSeqRecord
from gcgc.tests.fixtures import P53_HUMAN


class TestSeqRecordEncoder(unittest.TestCase):
    def setUp(self):
        self.sr = SeqRecord(Seq("ATCG"))
        self.alphabet = UnambiguousDnaAlphabet()

    def test_pad(self):
        padding = 10
        encoded = EncodedSeqRecord(self.alphabet, self.sr, padding_to=padding)
        self.assertEqual(len(encoded.pad), padding)

    def test_seq_record_encoder(self):
        expected_array = np.array(
            [
                [0, 0, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0],
            ],
            dtype=np.int,
        )
        encoded = EncodedSeqRecord(self.alphabet, self.sr)
        np.testing.assert_array_equal(encoded.one_hot_encode_sequence, expected_array)

    def test_yield_fasta_record(self):

        protein_alphabet = ExtendedProteinAlphabet()

        with open(P53_HUMAN) as f:
            for r in SeqIO.parse(f, format="fasta"):
                encoded = EncodedSeqRecord(protein_alphabet, r)
                self.assertEqual(
                    np.asarray(encoded.one_hot_encode_sequence).shape, (395, 29)
                )
