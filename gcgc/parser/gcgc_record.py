# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from dataclasses import dataclass
from pathlib import Path

from Bio.SeqRecord import SeqRecord


@dataclass
class GCGCRecord:
    path: Path
    seq_record: SeqRecord
