"""Unit test package for gcgc."""

import unittest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np

from gcgc.alphabet import UnambiguousDnaExtendedAlphabet
from gcgc.alphabet import IUPACProteinExtendedAlphabet
from gcgc.fasta_generator import FASTARecordGenerator
from gcgc.tests.fixtures import P53_HUMAN


class TestSeqRecordEncoder(unittest.TestCase):
    def setUp(self):
        self.sr = SeqRecord(Seq("ATCG"))
        self.alphabet = UnambiguousDnaExtendedAlphabet(5)

    def test_pad(self):
        alpha_len = 2
        alphabet = UnambiguousDnaExtendedAlphabet(alpha_len)
        encoded = alphabet.encode(self.sr)

        self.assertEqual(len(encoded), alpha_len)

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
        encoded_seq = self.alphabet.one_hot_encode_sequence(str(self.sr.seq))
        np.testing.assert_array_equal(encoded_seq, expected_array)

    def test_yield_fasta_record(self):

        protein_alphabet = IUPACProteinExtendedAlphabet(padding_to=50)
        fr = FASTARecordGenerator(P53_HUMAN, protein_alphabet)

        for encoded in fr:
            self.assertEqual(
                encoded["one_hot_sequence"].shape,
                (395, len(protein_alphabet.letters_and_tokens)),
            )
