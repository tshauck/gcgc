# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest
from pathlib import Path

import torch
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gcgc.fields import FileMetaDataField
from gcgc.parser.base import EncodedSeqLengthParser, SequenceParser


class TestParser(unittest.TestCase):
    def test_parser(self):
        vocab = [Path("ecoli"), Path("human")]

        f = FileMetaDataField.from_path_vocabulary("species", vocab)
        ff = [f]

        length_parser = EncodedSeqLengthParser(conform_to=10)

        sp = SequenceParser(encapsulate=True, seq_length_parser=length_parser, file_features=ff)

        dna = IUPAC.IUPACUnambiguousDNA()
        input_seq = SeqRecord(Seq("ATCG", alphabet=dna))

        test_values = [
            (input_seq, Path("ecoli"), torch.tensor(0)),
            (input_seq, Path("human"), torch.tensor(1)),
            (input_seq, Path("human"), torch.tensor(1)),
        ]

        for i, p, es in test_values:
            resp = sp.parse_record(i, p)
            self.assertEqual(resp["species"], es)
