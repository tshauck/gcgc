# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Objects and methods for dealing with PyTorch data."""

from pathlib import Path
from typing import Dict, Sequence

import torch
import torch.utils.data
from Bio import File, SeqIO

from gcgc.alphabet import ExtendedIUPACDNAEncoding
from gcgc.alphabet.base import EncodingAlphabet
from gcgc.ml.pytorch_utils.parser import TorchSequenceParser
from gcgc.parser.gcgc_record import GCGCRecord


class _SequenceIndexer(object):
    """A helper object that is used to index multiple files at ones."""

    def __init__(self):
        """Initialize the _SequenceIndexer object."""
        self._record_index = {}
        self._counter = -1

    def __call__(self, sid):
        """Return the record index if it's known, otherwise generate a new one."""
        try:
            return self._record_index[sid]
        except KeyError:
            self._counter = self._counter + 1
            self._record_index[sid] = self._counter
            return self._record_index[sid]


class GenomicDataset(torch.utils.data.Dataset):
    """GenomicDataset can be used to load sequence information into a format aminable to PyTorch."""

    def __init__(
        self, file_index: Dict[Path, File._IndexedSeqFileDict], parser: TorchSequenceParser
    ):
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
        file_format="fasta",
        alphabet: EncodingAlphabet = ExtendedIUPACDNAEncoding(),
    ) -> "GenomicDataset":
        """Initialize the GenomicDataset from a pathlib.Path sequence."""

        file_index = {}
        si = _SequenceIndexer()

        for f in sorted(path_sequence):
            file_index[f] = SeqIO.index(str(f), file_format, key_function=si, alphabet=alphabet)

        return cls(file_index, parser)

    def __len__(self) -> int:
        """Return the length of the dataset."""

        return sum(len(v) for v in self._file_index.values())

    def __getitem__(self, i: int):
        """Get the record from the index."""

        for k, v in self._file_index.items():
            try:
                r = GCGCRecord(path=k, seq_record=v[i])
                return self._parser.parse_record(r)
            except KeyError:
                pass

        raise RuntimeError(f"Exausted file index while looking for {i}.")
