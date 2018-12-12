# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest
from operator import itemgetter
from pathlib import Path

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gcgc.fields import AnnotationField, FileMetaDataField
from gcgc.parser.base import EncodedSeqLengthParser, SequenceParser
from gcgc.parser.gcgc_record import GCGCRecord


class TestParser(unittest.TestCase):
    def test_parser(self):
        vocab = [Path("ecoli"), Path("human")]

        f = FileMetaDataField.from_paths("species", vocab)
        ff = [f]

        annotations = {"annotation": "a"}
        af = AnnotationField.from_annotations(
            "annotation", [annotations], preprocess=itemgetter("annotation")
        )
        afs = [af]

        length_parser = EncodedSeqLengthParser(conform_to=10)

        sp = SequenceParser(
            encapsulate=True,
            seq_length_parser=length_parser,
            file_features=ff,
            annotation_features=afs,
        )

        dna = IUPAC.IUPACUnambiguousDNA()

        input_seq = SeqRecord(Seq("ATCG", alphabet=dna), annotations=annotations)

        test_values = [
            (input_seq, Path("ecoli"), 0, 0),
            (input_seq, Path("human"), 1, 0),
            (input_seq, Path("human"), 1, 0),
        ]

        for i, p, es, annotation_label in test_values:
            r = GCGCRecord(path=p, seq_record=i)
            resp = sp.parse_record(r)

            self.assertEqual(resp["species"], es)
            self.assertEqual(resp["annotation"], annotation_label)
