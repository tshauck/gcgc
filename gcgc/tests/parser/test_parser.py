# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from operator import itemgetter
from pathlib import Path

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gcgc.fields import AnnotationField, DescriptionField, FileMetaDataField
from gcgc.alphabet import IUPACUnambiguousDNAEncoding
from gcgc.parser.base import EncodedSeqLengthParser, SequenceParser
from gcgc.parser.gcgc_record import GCGCRecord


def test_parser():
    def preprocess(d: str) -> str:
        return d.split("\t")[-1]

    d = DescriptionField.from_descriptions("desc_field", ["Test\tA", "Test\tB"], preprocess)

    kmer_step_size = 2

    vocab = [Path("ecoli"), Path("human")]
    f = FileMetaDataField.from_paths("species", vocab)

    annotations = {"annotation": "a"}
    af = AnnotationField.from_annotations(
        "annotation", [annotations], preprocess=itemgetter("annotation")
    )

    length_parser = EncodedSeqLengthParser(conform_to=5)

    dna = IUPACUnambiguousDNAEncoding(kmer_size=kmer_step_size)

    sp = SequenceParser(
        encapsulate=True,
        seq_length_parser=length_parser,
        file_features=[f],
        annotation_features=[af],
        description_features=[d],
        kmer_step_size=kmer_step_size,
    )

    test_values = [
        ("ATCGATCG", Path("ecoli"), 0, 0, 0, [1, 9, 15, 0, 0]),
        ("ACT", Path("human"), 1, 0, 0, [1, 10, 2, 0, 0]),
        ("ATCGATGG", Path("human"), 1, 0, 0, [1, 9, 15, 0, 0]),
    ]

    for s, p, es, annotation_label, description_label, expected_seq in test_values:

        i = SeqRecord(Seq(s, alphabet=dna), annotations=annotations, description="Test\tA")
        r = GCGCRecord(path=p, seq_record=i)
        resp = sp.parse_record(r, parsed_seq_len=5)

        assert resp["species"] == es
        assert resp["annotation"] == annotation_label
        assert resp["desc_field"] == description_label
        assert resp["seq_tensor"] == expected_seq
