# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A GCGC Record object."""

from pathlib import Path

from Bio.SeqRecord import SeqRecord


class GCGCRecord:
    """A class for holding the Path of the file and a SeqRecord."""

    def __init__(self, path: Path, seq_record: SeqRecord) -> None:
        self.path = path
        self.seq_record = seq_record
