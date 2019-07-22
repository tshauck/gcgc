# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from operator import itemgetter
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gcgc.alphabet import IUPACUnambiguousDNAEncoding
from gcgc.fields import AnnotationField
from gcgc.fields import DescriptionField
from gcgc.fields import FileMetaDataField
from gcgc.parser.base import EncodedSeqLengthParser
from gcgc.parser.base import SequenceParser
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

    length_parser = EncodedSeqLengthParser(conform_to=10)

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
        ("ATCGATCGATCGATCGATCG", Path("ecoli"), 0, 0, 0, [1, 9, 15, 9, 15], 5),
        ("ACT", Path("human"), 1, 0, 0, [1, 10, 2, 0, 0], 3),
        ("ATCGATGG", Path("human"), 1, 0, 0, [1, 9, 15, 9, 3], 5),
    ]

    for s, p, es, annotation_label, description_label, expected_seq, expected_len in test_values:

        i = SeqRecord(Seq(s, alphabet=dna), annotations=annotations, description="Test\tA")
        r = GCGCRecord(path=p, seq_record=i)
        resp = sp.parse_record(r, parsed_seq_len=5)

        assert resp["species"] == es
        assert resp["annotation"] == annotation_label
        assert resp["desc_field"] == description_label
        assert resp["seq_tensor"] == expected_seq
        assert resp["seq_len"] == expected_len


def test_masked_parser():

    sp = SequenceParser(masked_probability=0.5)
    test_values = [("ATCGATCGATCGATCGATCG", []), ("ACTATCAATT", []), ("ATCGATGG", [])]
    dna = IUPACUnambiguousDNAEncoding(masked=True, kmer_size=2)

    for seq, out in test_values:
        i = SeqRecord(Seq(seq, alphabet=dna))
        r = GCGCRecord(seq_record=i, path=None)
        record = sp.parse_record(r)

        assert "seq_tensor_masked" in record
        assert (record["seq_tensor_masked"] == dna.encoded_mask).sum() > 0
