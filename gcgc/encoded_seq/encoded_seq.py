# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Contains the EncodedSeq object."""

from typing import Iterable, Sequence, Union

from Bio.Seq import Seq
import numpy as np

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

        if seq.alphabet and isinstance(seq.alphabet, EncodingAlphabet):
            return cls(str(seq), seq.alphabet)

        gcgc_alphabet = biopython_alphabet_to_gcgc_alphabet(seq.alphabet, *args, **kwargs)
        return cls(str(seq), gcgc_alphabet)

    @property
    def has_start_token(self) -> bool:
        """Return True if the sequence contains the start token."""
        return self.alphabet.START in self

    @property
    def has_end_token(self) -> bool:
        """Return True if the sequence contains the end token."""
        return self.alphabet.END in self

    def encapsulate(self, aware: bool = True) -> "EncodedSeq":
        """Return a sequence with the alphabet start and end token applied to the start and end."""

        if aware and self.has_start_token:
            start_char = ""
        else:
            start_char = self.alphabet.START

        if aware and self.has_end_token:
            end_char = ""
        else:
            end_char = self.alphabet.END

        return EncodedSeq.from_seq(start_char + self + end_char, self.alphabet)

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
        if seq_len < conform_to:
            return self.pad(pad_to=conform_to)

        return self[:conform_to]

    def shift(self, offset: int) -> "EncodedSeq":
        """Shift the by the offset and return a new sequence."""

        if offset == 0:
            return self

        if offset > 0:
            padding_characters = (offset - 1) * self.alphabet.PADDING

            if self.has_start_token:
                start_char = self.alphabet.PADDING
            else:
                start_char = self.alphabet.START

            starting_characters = padding_characters + start_char
            return starting_characters + self[: (-1 * offset)]

        if offset < 0:
            offset_abs = abs(offset)
            padding_characters = (offset_abs - 1) * self.alphabet.PADDING

            if self.has_end_token:
                end_char = self.alphabet.PADDING
            else:
                end_char = self.alphabet.END

            ending_characters = end_char + padding_characters
            return self[offset_abs:] + ending_characters

        raise ValueError(f"Unsure how to handle {offset}.")

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

    @classmethod
    def from_integer_encoded_seq(
        cls, integer_encoded_seq: Iterable[int], alphabet: EncodingAlphabet
    ) -> "EncodedSeq":
        """Given the encoded seq and alphabet, return the EncodedSeq."""

        seq_str = alphabet.decode_tokens(integer_encoded_seq)
        return cls(seq_str, alphabet)

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
