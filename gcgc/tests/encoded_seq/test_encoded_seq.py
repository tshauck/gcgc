# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from numpy.testing import assert_array_equal
import pytest

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

    def test_one_hot_encoding(self):

        expected = [
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
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


def test_decode_tokens():
    alphabet = ExtendedIUPACDNAEncoding()

    new_seq = EncodedSeq.from_integer_encoded_seq([1, 4, 4, 2], alphabet)
    expected_seq = EncodedSeq(">AA<", alphabet)

    assert new_seq == expected_seq


@pytest.mark.parametrize(
    "actual_seq,expected_seq,offset",
    [
        (
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()),
            EncodedSeq(">ATC", ExtendedIUPACDNAEncoding()),
            1,
        ),
        (
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()).encapsulate(),
            EncodedSeq("|>ATCG", ExtendedIUPACDNAEncoding()),
            1,
        ),
        (
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()),
            EncodedSeq("|>AT", ExtendedIUPACDNAEncoding()),
            2,
        ),
        (
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()),
            EncodedSeq("TCG<", ExtendedIUPACDNAEncoding()),
            -1,
        ),
        (
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()).encapsulate(),
            EncodedSeq("ATCG<|", ExtendedIUPACDNAEncoding()),
            -1,
        ),
        (
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()),
            EncodedSeq("CG<|", ExtendedIUPACDNAEncoding()),
            -2,
        ),
        (
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()),
            EncodedSeq("ATCG", ExtendedIUPACDNAEncoding()),
            0,
        ),
    ],
)
def test_seq_shift(actual_seq, expected_seq, offset):
    shifted = actual_seq.shift(offset)
    assert shifted == expected_seq


@pytest.mark.parametrize(
    "seq,kmer_size,integer_encoding", [(">ATCG<||", 2, [1, 13, 22, 27, 2, 0, 0])]
)
def test_kmer_size_encoding(seq, kmer_size, integer_encoding):
    encoded_seq = EncodedSeq(seq, ExtendedIUPACDNAEncoding(kmer_size=kmer_size))
    assert encoded_seq.integer_encoded == integer_encoding
