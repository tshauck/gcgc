# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest
from typing import NamedTuple, Dict, List, Tuple

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from numpy.testing import assert_array_equal

from gcgc.alphabet.iupac import ExtendedIUPACDNAEncoding
from gcgc.encoded_seq import EncodedSeq
from gcgc.exceptions import GCGCAlphabetException


class TestEncodedSeq(unittest.TestCase):
    def test_raises_for_bad_alphabet(self):
        with self.assertRaises(GCGCAlphabetException):
            EncodedSeq("ATCG", 1)

    def test_pad_to(self):

        pad_to = 10
        es = EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()).pad(pad_to)

        self.assertEqual(len(es), pad_to)

    def test_pad_to_over(self):

        pad_to = 2
        es = EncodedSeq("ATCG", ExtendedIUPACDNAEncoding())
        new_es = es.pad(pad_to)

        self.assertEqual(len(new_es), len(es))

    def test_encapsulate(self):

        es = EncodedSeq("ATCG", ExtendedIUPACDNAEncoding())
        new_es = es.encapsulate()
        self.assertEqual(new_es[0], new_es.alphabet.START)
        self.assertEqual(new_es[-1], new_es.alphabet.END)

    def test_conform(self):

        length = 5

        es = EncodedSeq("A", ExtendedIUPACDNAEncoding())
        new_es = es.encapsulate().conform(length)
        self.assertEqual(len(new_es), length)
        self.assertIsInstance(new_es, EncodedSeq)

        es = EncodedSeq("ATCGGCG", ExtendedIUPACDNAEncoding())
        new_es = es.encapsulate().conform(length)
        self.assertEqual(len(new_es), length)
        self.assertIsInstance(new_es, EncodedSeq)

        es = EncodedSeq("ATC", ExtendedIUPACDNAEncoding())
        new_es = es.encapsulate().conform(length)
        self.assertEqual(len(new_es), length)
        self.assertIsInstance(new_es, EncodedSeq)

    def test_kmer_rollout(self):
        class KMerTestSet(NamedTuple):
            name: str
            encoded_seq: EncodedSeq
            rollout_options: Dict[str, int]
            expected_seqs: List[Tuple[str, str, str]]

        test_table = [
            KMerTestSet(
                name="Test length of kmer",
                encoded_seq=EncodedSeq("ATCGATCG", ExtendedIUPACDNAEncoding()),
                rollout_options={"kmer_length": 4},
                expected_seqs=[
                    ("ATCG", "", "A"),
                    ("TCGA", "", "T"),
                    ("CGAT", "", "C"),
                    ("GATC", "", "G"),
                ],
            ),
            KMerTestSet(
                name="Test length and prior length.",
                encoded_seq=EncodedSeq("ATCGATCG", ExtendedIUPACDNAEncoding()),
                rollout_options={"kmer_length": 4, "prior_length": 2},
                expected_seqs=[("CGAT", "AT", "C"), ("GATC", "TC", "G")],
            ),
            KMerTestSet(
                name="Test Window",
                encoded_seq=EncodedSeq("ATCGATCGATCG", ExtendedIUPACDNAEncoding()),
                rollout_options={"kmer_length": 4, "prior_length": 2, "window": 2},
                expected_seqs=[("CGAT", "AT", "C"), ("ATCG", "CG", "A"), ("CGAT", "AT", "C")],
            ),
        ]

        for test_set in test_table:
            rollout_iters = test_set.encoded_seq.rollout_kmers(**test_set.rollout_options)
            for i, rk in enumerate(rollout_iters):

                expected_kmer, expected_prior, expected_next = test_set.expected_seqs[i]

                self.assertEqual(str(expected_kmer), rk.kmer)
                self.assertEqual(str(expected_prior), rk.prior_kmer)
                self.assertEqual(str(expected_next), rk.next_kmer)

    def test_one_hot_encoding(self):

        expected = [
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]

        es = EncodedSeq("ATCG", ExtendedIUPACDNAEncoding())
        assert_array_equal(es.one_hot_encoded, expected)

    def test_from_seq_bad_alphabet(self):
        seq = Seq("ATCG", None)

        with self.assertRaises(GCGCAlphabetException):
            EncodedSeq.from_seq(seq)

    def test_from_seq(self):

        alphabet = IUPAC.IUPACUnambiguousDNA()
        seq = Seq("ATCG", alphabet)
        encoded_seq = EncodedSeq.from_seq(seq)

        self.assertEqual(seq, encoded_seq)
