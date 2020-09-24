# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Module holding pytorch code."""

from pathlib import Path
from typing import Dict
from typing import Sequence

from gcgc.tokenizer.kmer_tokenzier import KmerTokenizer

try:
    import torch
    import torch.utils.data

    from Bio import File
    from Bio import SeqIO

except ImportError as exp:
    # pylint: disable=invalid-name
    needed = "torch, biopython"
    raise ImportError(f"Missing one or more libraries: {needed}. Please install: {exp}") from exp


# pylint: disable=too-few-public-methods
class _SequenceIndexer:
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

    def __init__(self, file_index: Dict[Path, File._IndexedSeqFileDict], tokenizer: KmerTokenizer):
        """Initialize the GenomicDataset object."""
        self._file_index = file_index
        self._tokenizer = tokenizer

        super().__init__()

    @classmethod
    def from_path(
        cls, path: Path, tokenizer: KmerTokenizer, file_format: str = "fasta",
    ) -> "GenomicDataset":
        """Init from a single file. This is a convenience method that delegates to from_paths."""
        return cls.from_paths([path], tokenizer, file_format)

    @classmethod
    def from_paths(
        cls, path_sequence: Sequence[Path], tokenizer: KmerTokenizer, file_format="fasta",
    ) -> "GenomicDataset":
        """Initialize the GenomicDataset from a pathlib.Path sequence."""
        file_index = {}
        seq_indexer = _SequenceIndexer()

        for file_path in sorted(path_sequence):
            file_index[file_path] = SeqIO.index(
                str(file_path), file_format, key_function=seq_indexer
            )

        return cls(file_index, tokenizer)

    def __len__(self) -> int:
        """Return the length of the dataset."""
        return sum(len(v) for v in self._file_index.values())

    def __getitem__(self, i: int):
        """Get the record from the index."""
        for value in self._file_index.values():
            try:
                seq_record = value[i]
            except KeyError:
                continue

            tokenized = self._tokenizer.encode(str(seq_record.seq))
            # pylint: disable=not-callable, no-member
            return torch.tensor(tokenized, dtype=torch.long)

        # pylint: disable=useless-else-on-loop
        else:
            raise RuntimeError(f"Exausted file index while looking for {i}.")