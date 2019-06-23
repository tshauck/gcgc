# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Objects and methods for dealing with PyTorch data."""

from pathlib import Path
from typing import Sequence

from Bio import File
from Bio import SeqIO
import torch
import torch.utils.data

from gcgc.alphabet import ExtendedIUPACDNAEncoding
from gcgc.alphabet.base import EncodingAlphabet
from gcgc.ml.pytorch_utils.parser import TorchSequenceParser
from gcgc.parser.gcgc_record import GCGCRecord


class GenomicDataset(torch.utils.data.Dataset):
    """GenomicDataset can be used to load sequence information into a format aminable to PyTorch."""

    def __init__(self, file_index: File._SQLiteManySeqFilesDict, parser: TorchSequenceParser):
        """Initialize the GenomicDataset object."""

        self._file_index = file_index
        self._parser = parser

        super().__init__()

    @classmethod
    def from_path(
        cls,
        path: Path,
        parser: TorchSequenceParser,
        file_format: str = "fasta",
        alphabet: EncodingAlphabet = ExtendedIUPACDNAEncoding(),
    ) -> "GenomicDataset":
        """Init from a single file. This is a convience method that delegates to from_paths."""

        return cls.from_paths([path], parser, file_format, alphabet)

    @classmethod
    def from_paths(
        cls,
        path_sequence: Sequence[Path],
        parser: TorchSequenceParser,
        file_format: str = "fasta",
        alphabet: EncodingAlphabet = ExtendedIUPACDNAEncoding(),
        index_db: str = ":memory:",
        **kwargs,
    ) -> "GenomicDataset":
        """Initialize the GenomicDataset from a pathlib.Path sequence."""

        file_index = SeqIO.index_db(
            index_db, [str(p) for p in path_sequence], file_format, alphabet=alphabet, **kwargs
        )
        return cls(file_index, parser)

    def __len__(self) -> int:
        """Return the length of the dataset."""

        return len(self._file_index)

    def __getitem__(self, i: int):
        """Get the record from the index."""

        qry = "SELECT key, file_number FROM offset_data LIMIT 1 OFFSET ?;"
        key, file_number = self._file_index._con.execute(qry, (i,)).fetchone()
        file_name = Path(self._file_index._filenames[file_number])

        r = GCGCRecord(path=file_name, seq_record=self._file_index[key])
        return self._parser.parse_record(r)
