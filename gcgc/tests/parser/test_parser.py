# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from operator import itemgetter
from pathlib import Path

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gcgc.fields import AnnotationField, DescriptionField, FileMetaDataField
from gcgc.parser.base import EncodedSeqLengthParser, SequenceParser
from gcgc.parser.gcgc_record import GCGCRecord


def test_parser():
    def preprocess(d: str) -> str:
        return d.split("\t")[-1]

    d = DescriptionField.from_descriptions("desc_field", ["Test\tA", "Test\tB"], preprocess)

    vocab = [Path("ecoli"), Path("human")]
    f = FileMetaDataField.from_paths("species", vocab)

    annotations = {"annotation": "a"}
    af = AnnotationField.from_annotations(
        "annotation", [annotations], preprocess=itemgetter("annotation")
    )

    length_parser = EncodedSeqLengthParser(conform_to=10)

    sp = SequenceParser(
        encapsulate=True,
        seq_length_parser=length_parser,
        file_features=[f],
        annotation_features=[af],
        description_features=[d],
    )

    dna = IUPAC.IUPACUnambiguousDNA()

    input_seq = SeqRecord(Seq("ATCG", alphabet=dna), annotations=annotations, description="Test\tA")

    test_values = [
        (input_seq, Path("ecoli"), 0, 0, 0),
        (input_seq, Path("human"), 1, 0, 0),
        (input_seq, Path("human"), 1, 0, 0),
    ]

    for i, p, es, annotation_label, description_label in test_values:
        r = GCGCRecord(path=p, seq_record=i)
        resp = sp.parse_record(r)

        assert resp["species"] == es
        assert resp["annotation"] == annotation_label
        assert resp["desc_field"] == description_label
