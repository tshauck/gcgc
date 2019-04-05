# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from pathlib import Path

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import torch

from gcgc.fields import FileMetaDataField
from gcgc.ml.pytorch_utils.parser import TorchSequenceParser
from gcgc.parser.base import EncodedSeqLengthParser
from gcgc.parser.gcgc_record import GCGCRecord


def test_parser():
    vocab = [Path("ecoli"), Path("human")]

    f = FileMetaDataField.from_paths("species", vocab)
    ff = [f]

    length_parser = EncodedSeqLengthParser(conform_to=10)

    sp = TorchSequenceParser(
        encapsulate=True, seq_length_parser=length_parser, file_features=ff, sequence_offset=-1
    )

    dna = IUPAC.IUPACUnambiguousDNA()
    input_seq = SeqRecord(Seq("ATCG", alphabet=dna))

    test_values = [
        (input_seq, Path("ecoli"), torch.tensor(0)),
        (input_seq, Path("human"), torch.tensor(1)),
        (input_seq, Path("human"), torch.tensor(1)),
    ]

    for i, p, es in test_values:
        r = GCGCRecord(path=p, seq_record=i)
        resp = sp.parse_record(r)
        assert resp["species"] == es
