# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from pathlib import Path
from typing import Dict, Sequence

import torch
import torch.utils.data
from Bio import File, SeqIO

from gcgc.alphabet import ExtendedIUPACDNAEncoding
from gcgc.alphabet.base import EncodingAlphabet
from gcgc.parser.gcgc_record import GCGCRecord
from gcgc.third_party.pytorch_utils.parser import TorchSequenceParser


class _SequenceIndexer(object):
    """
    This is a helper object that is used to mediate between the indexes in the underlying files
    and providing a consistent interface for __getitem__ in the GenomicDataset object.
    """

    def __init__(self):
        self._record_index = {}
        self._counter = -1

    def __call__(self, sid):
        try:
            return self._record_index[sid]
        except KeyError:
            self._counter = self._counter + 1
            self._record_index[sid] = self._counter
            return self._record_index[sid]


class GenomicDataset(torch.utils.data.Dataset):
    """
    GenomicDataset can be used to load sequence information into a format aminable to PyTorch.
    """

    def __init__(
        self,
        file_index: Dict[Path, File._IndexedSeqFileDict],
        parser: TorchSequenceParser,
        file_format: str = "fasta",
    ):
        """
        Initialize the GenomicDataset object.

        It is recommended that you do _not_ itstantiate this object directly, and instead you should
        use one of the @classmethods.

        Args:
            file_index: A dictionary mapping the Path of the file to its index.
            parser: The parser object which specifies how the input sequence should be parsed.
            file_format: The format of the file to parse. This should be one understood by
                Biopython.
        """

        self._file_index = file_index
        self._parser = parser
        self._file_format = file_format

        super().__init__()

    @classmethod
    def from_path(
        cls,
        path: Path,
        parser: TorchSequenceParser,
        file_format: str = "fasta",
        alphabet: EncodingAlphabet = ExtendedIUPACDNAEncoding(),
    ) -> "GenomicDataset":
        """
        Init from a single file. This is a convience method that delegates to
        from_paths.
        """
        return cls.init_from_path_generator([path], parser, file_format, alphabet)

    @classmethod
    def from_paths(
        cls,
        path_sequence: Sequence[Path],
        parser: TorchSequenceParser,
        file_format="fasta",
        alphabet: EncodingAlphabet = ExtendedIUPACDNAEncoding(),
    ) -> "GenomicDataset":
        """
        Initialize the GenomicDataset from a pathlib.Path sequence. For example, the
        pathlib.Path.glob returns a Path generator.
        """

        file_index = {}
        si = _SequenceIndexer()

        for f in sorted(path_sequence):
            file_index[f] = SeqIO.index(str(f), file_format, key_function=si, alphabet=alphabet)

        return cls(file_index, parser=parser, file_format=file_format)

    def __len__(self) -> int:
        """
        Return the length of the dataset.
        """
        return sum(len(v) for v in self._file_index.values())

    def __getitem__(self, i: int):
        """
        Get the record from the index.
        """

        for k, v in self._file_index.items():
            try:
                r = GCGCRecord(path=k, seq_record=v[i])
                return self._parser.parse_record(r)
            except KeyError:
                pass

        else:
            raise RuntimeError(f"Exausted file index while looking for {i}.")
