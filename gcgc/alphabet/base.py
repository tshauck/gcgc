# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Holds the base EncodingAlphabet."""

from abc import ABC
from typing import Sequence

from gcgc.exceptions import GCGCAlphabetLetterEncodingException


class EncodingAlphabet(ABC):
    """The Encoding Alphabet is meant to be a baseclass for other alphabets."""

    START: str = ">"
    END: str = "<"
    PADDING: str = "|"

    def __init__(self):
        """Create the EncodingAlphabet object."""

        self.letters_and_tokens = self.letters + self.START + self.END + self.PADDING

        self.encoding_index = {letter: idx for idx, letter in enumerate(self.letters_and_tokens)}
        self.decoding_index = {idx: letter for letter, idx in self.encoding_index.items()}

    def __len__(self) -> int:
        """Get the lenght of the Alphabet."""
        return len(self.letters_and_tokens)

    def encode_token(self, token: str) -> int:
        """Given a particular token, return the integer representation."""

        return self.encoding_index[token]

    def decode_token(self, int_token: int) -> str:
        """Decode a token. This is the inverse of encode_token."""

        return self.decoding_index[int_token]

    def integer_encode(self, seq: str) -> Sequence[int]:
        """Integer encode the sequence."""

        try:
            return [self.encoding_index[s] for s in seq]
        except KeyError:
            raise GCGCAlphabetLetterEncodingException(f"{seq} not in {self.encoding_index}")

    def integer_decode(self, int_seq: Sequence[int]) -> str:
        """Given a sequence of integers, convert it to a string."""

        return "".join(self.decode_token(s) for s in int_seq)
