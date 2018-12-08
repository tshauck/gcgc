# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from typing import Sequence

import numpy as np
from Bio.Seq import Seq

from gcgc.alphabet.base import EncodingAlphabet
from gcgc.alphabet.utils import biopython_alphabet_to_gcgc_alphabet
from gcgc.exceptions import GCGCAlphabetException


class EncodedSeq(Seq):
    def __init__(self, data, alphabet) -> None:
        if not isinstance(alphabet, EncodingAlphabet):
            raise GCGCAlphabetException(f"Cannot use alphabet of type {type(alphabet)}.")

        super().__init__(data, alphabet)

    @classmethod
    def from_seq(cls, seq: Seq):
        """
        Instantiate an EncodedSeq object from a Seq object.
        """

        gcgc_alphabet = biopython_alphabet_to_gcgc_alphabet(seq.alphabet)
        return cls(str(seq), gcgc_alphabet)

    def encapsulate(self) -> "EncodedSeq":
        return self.alphabet.START + self + self.alphabet.END

    def pad(self, pad_to: int = 50) -> "EncodedSeq":
        """
        Pad a sequence up to `pad_to` characters.
        """

        seq_len = len(self)

        if seq_len < pad_to:
            n_extra_chars = pad_to - seq_len
            extra_chars = self.alphabet.PADDING * n_extra_chars
        else:
            extra_chars = ""

        return self + extra_chars

    def conform(self, conform_to: int = 50) -> "EncodedSeq":
        seq_len = len(self)

        if seq_len == conform_to:
            return self
        elif seq_len < conform_to:
            return self.pad(pad_to=conform_to)
        else:
            return self[:conform_to]

    @property
    def integer_encoded(self):
        return self.alphabet.integer_encode(self)

    @property
    def one_hot_encoded(self) -> Sequence[Sequence[int]]:
        """
        Encodes D x N where D is the size of the alphabet and N is the padding.

        Returns:
            A one hot encoded matrix representing the sequence.
        """

        encoded_sequence = self.alphabet.integer_encode(self)
        encoded_len = len(encoded_sequence)
        letters_len = len(self.alphabet.letters_and_tokens)

        one_hot_seq = np.zeros((encoded_len, letters_len), dtype=np.int)
        one_hot_seq[np.arange(encoded_len), encoded_sequence] = 1

        return one_hot_seq.tolist()

    def __add__(self, other) -> "EncodedSeq":
        """
        Add two enccoded sequences together.
        """
        added_seq = super().__add__(other)
        return self.from_seq(added_seq)

    def __getitem__(self, index) -> "EncodedSeq":
        got_item = super().__getitem__(index)
        if isinstance(index, int):
            return got_item
        else:
            return self.from_seq(got_item)
