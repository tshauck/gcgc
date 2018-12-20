# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A GCGC Record object."""

from pathlib import Path

from Bio.SeqRecord import SeqRecord

from gcgc.encoded_seq import EncodedSeq


class GCGCRecord(object):
    """A class for holding the Path of the file and a SeqRecord."""

    def __init__(self, path: Path, seq_record: SeqRecord) -> None:
        self.path = path
        self.seq_record = seq_record

    @property
    def encoded_seq(self) -> EncodedSeq:
        return EncodedSeq.from_seq(self.seq_record.seq, self.seq_record.seq.alphabet)
