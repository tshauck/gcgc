# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Module holding third party code."""

import json
import os.path
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Sequence
from typing import Tuple

import torch
import torch.utils.data

from Bio import File
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from transformers.tokenization_utils import PreTrainedTokenizer

from gcgc.tokenizer.kmer_tokenzier import KmerTokenizer


class GCGCTransformersTokenizer(PreTrainedTokenizer):
    """A GCGC Tokenizer that is compatible with the transformers library."""

    VOCABULARY_FILENAME = "gcgc-vocab.json"

    def __init__(self, kmer_tokenizer: KmerTokenizer, **kwargs: Dict[str, Any]):
        """Init the GCGCTransformersTokenizer object."""
        self.kmer_tokenizer = kmer_tokenizer

        super().__init__(
            bos_token=self.kmer_tokenizer.bos_token,
            eos_token=self.kmer_tokenizer.eos_token,
            unk_token=self.kmer_tokenizer.unk_token,
            pad_token=self.kmer_tokenizer.pad_token,
            mask_token=self.kmer_tokenizer.mask_token,
            **kwargs,
        )

    def get_vocab(self):
        """Return the vocabulary for this transformer."""
        return self.kmer_tokenizer.vocab.stoi

    @property
    def vocab_size(self):
        """Return the size of the vocabulary."""
        return len(self.kmer_tokenizer.vocab)

    def _convert_id_to_token(self, index: int) -> str:
        return self.kmer_tokenizer.vocab.itos[index]

    def _convert_token_to_id(self, token: str) -> int:
        return self.kmer_tokenizer.vocab.stoi[token]

    def _tokenize(self, text: str, **kwargs):
        return self.kmer_tokenizer.tokenize(text)

    def convert_tokens_to_string(self, tokens: List[str]) -> str:
        """Convert the tokens into a string."""
        return "".join(tokens)

    def save_vocabulary(self, save_directory) -> Tuple[str]:
        """Save the vocabulary string to integer map in the save_directory."""
        vocab_location = os.path.join(save_directory, self.VOCABULARY_FILENAME)

        with open(vocab_location, "w") as outf:
            json.dump(self.kmer_tokenizer.vocab.stoi, outf)

        return (vocab_location,)


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

    def __init__(
        self, file_index: Dict[Path, File._IndexedSeqFileDict], tokenizer: GCGCTransformersTokenizer
    ):
        """Initialize the GenomicDataset object."""
        self._file_index = file_index
        self._tokenizer = tokenizer

        super().__init__()

    @classmethod
    def from_path(
        cls,
        path: Path,
        tokenizer: GCGCTransformersTokenizer,
        file_format: str = "fasta",
        alphabet=IUPAC.ExtendedIUPACProtein(),
    ) -> "GenomicDataset":
        """Init from a single file. This is a convenience method that delegates to from_paths."""
        return cls.from_paths([path], tokenizer, file_format, alphabet)

    @classmethod
    def from_paths(
        cls,
        path_sequence: Sequence[Path],
        tokenizer: GCGCTransformersTokenizer,
        file_format="fasta",
        alphabet=IUPAC.ExtendedIUPACProtein(),
    ) -> "GenomicDataset":
        """Initialize the GenomicDataset from a pathlib.Path sequence."""
        file_index = {}
        seq_indexer = _SequenceIndexer()

        for file_path in sorted(path_sequence):
            file_index[file_path] = SeqIO.index(
                str(file_path), file_format, key_function=seq_indexer, alphabet=alphabet
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
