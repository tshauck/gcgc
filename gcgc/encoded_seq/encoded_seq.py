# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Contains the EncodedSeq object."""

from typing import Sequence, Union

import numpy as np
from Bio.Seq import Seq

from gcgc.alphabet.base import EncodingAlphabet
from gcgc.alphabet.utils import biopython_alphabet_to_gcgc_alphabet
from gcgc.exceptions import GCGCAlphabetException


class EncodedSeq(Seq):
    """A wrapper around Bio.Seq.Seq."""

    def __init__(self, data, alphabet) -> None:
        """Initialize the EncodedSeq object with sequence data and an alphabet."""

        if not isinstance(alphabet, EncodingAlphabet):
            raise GCGCAlphabetException(f"Cannot use alphabet of type {type(alphabet)}.")

        super().__init__(data, alphabet)

    @classmethod
    def from_seq(cls, seq: Seq, gcgc_alphabet=None, *args, **kwargs) -> "EncodedSeq":
        """Instantiate an EncodedSeq object from a Seq object."""

        if gcgc_alphabet and isinstance(gcgc_alphabet, EncodingAlphabet):
            return cls(str(seq), gcgc_alphabet)

        gcgc_alphabet = biopython_alphabet_to_gcgc_alphabet(seq.alphabet, *args, **kwargs)
        return cls(str(seq), gcgc_alphabet)

    def encapsulate(self) -> "EncodedSeq":
        """Return a sequence with the alphabet start and end token applied to the start and end."""
        return EncodedSeq.from_seq(self.alphabet.START + self + self.alphabet.END, self.alphabet)

    def pad(self, pad_to: int = 50) -> "EncodedSeq":
        """Pad a sequence up to `pad_to` characters."""

        seq_len = len(self)

        if seq_len < pad_to:
            n_extra_chars = pad_to - seq_len
            extra_chars = self.alphabet.PADDING * n_extra_chars
        else:
            extra_chars = ""

        return self + extra_chars

    def conform(self, conform_to: int = 50) -> "EncodedSeq":
        """Return a exactly equal to the conform_to value."""

        seq_len = len(self)

        if seq_len == conform_to:
            return self
        elif seq_len < conform_to:
            return self.pad(pad_to=conform_to)
        else:
            return self[:conform_to]

    @property
    def integer_encoded(self):
        """Return the underlying sequence in its integer representation."""
        return self.alphabet.integer_encode(self)

    @property
    def one_hot_encoded(self) -> Sequence[Sequence[int]]:
        """Encode into D x N matrix where D is the size of the alphabet and N is the padding."""

        encoded_sequence = self.alphabet.integer_encode(self)
        encoded_len = len(encoded_sequence)
        letters_len = len(self.alphabet.letters_and_tokens)

        one_hot_seq = np.zeros((encoded_len, letters_len), dtype=np.int)
        one_hot_seq[np.arange(encoded_len), encoded_sequence] = 1

        return one_hot_seq.tolist()

    def __add__(self, other: "EncodedSeq") -> "EncodedSeq":
        """Add two enccoded sequences together."""

        added_seq = super().__add__(other)
        return self.from_seq(added_seq, self.alphabet)

    def __getitem__(self, index: Union[int, slice]) -> "EncodedSeq":
        """Given the input index for the entire datset, return the associated EncodedSeq."""

        got_item = super().__getitem__(index)
        if isinstance(index, int):
            return got_item
        else:
            return self.from_seq(got_item)
