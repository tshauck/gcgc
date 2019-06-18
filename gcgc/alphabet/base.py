# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Holds the base EncodingAlphabet."""

import itertools as it
from typing import Iterable
from typing import Sequence

from gcgc.exceptions import GCGCAlphabetLetterEncodingException


class EncodingAlphabet:
    """The Encoding Alphabet is meant to be a baseclass for other alphabets."""

    START: str = ">"
    END: str = "<"
    PADDING: str = "|"

    # Convince linting that EncodingAlphabet will have a letters attribute.
    letters: str

    def __init__(self, kmer_size: int = 1):
        """Create the EncodingAlphabet object."""

        self.letters_and_tokens = self.START + self.END + self.PADDING + self.letters
        self.kmer_size = kmer_size

        self.kmers = ["".join(kmer) for kmer in it.product(self.letters, repeat=self.kmer_size)]
        self.kmers_and_tokens = list(self.START) + list(self.END) + list(self.PADDING) + self.kmers

        self.encoding_index = {letter: idx for idx, letter in enumerate(self.kmers_and_tokens)}
        self.decoding_index = {idx: letter for letter, idx in self.encoding_index.items()}

    @property
    def encoded_padding(self):
        """Get the integer for the padding character."""
        return self.encode_token(self.PADDING)

    @property
    def encoded_start(self):
        """Get the integer for the start character."""
        return self.encode_token(self.START)

    @property
    def encoded_end(self):
        """Get the integer for the end character."""
        return self.encode_token(self.END)

    def __len__(self) -> int:
        """Get the lenght of the Alphabet."""
        return len(self.letters_and_tokens)

    def encode_token(self, token: str) -> int:
        """Given a particular token, return the integer representation."""
        return self.encoding_index[token]

    def decode_token(self, int_token: int) -> str:
        """Decode a token. This is the inverse of encode_token."""
        return self.decoding_index[int_token]

    def decode_tokens(self, int_tokens: Iterable[int]) -> str:
        """Decode an iterable of integer tokens into a single string."""
        return "".join(self.decode_token(t) for t in int_tokens)

    def _kmer_one(self, seq):
        try:
            encoded = []
            seq_len = len(seq)

            for i in range(0, seq_len):
                kmer = seq[i]
                encoded.append(self.encoding_index[kmer])
            return encoded

        except KeyError:
            raise GCGCAlphabetLetterEncodingException(f"{kmer} not in {self.encoding_index}")

    def _kmer_n(self, seq: str) -> Sequence[int]:
        try:
            encoded = []

            seq_len = len(seq)
            iterations = seq_len - self.kmer_size + 1

            for i in range(0, iterations):
                kmer = seq[i : i + self.kmer_size]
                encoded.append(self.encoding_index[kmer])
            return encoded

        except KeyError:
            raise GCGCAlphabetLetterEncodingException(f"{kmer} not in {self.encoding_index}")

    def integer_encode(self, seq: str) -> Sequence[int]:
        """Integer encode the sequence."""

        stripped_seq = "".join(s for s in seq if s not in {self.START, self.END, self.PADDING})
        seq_len = len(stripped_seq)

        if seq_len < self.kmer_size:
            raise ValueError(
                f"seq length {seq_len} cannot be less than the kmer size {self.kmer_size}"
            )

        if self.kmer_size == 1:
            encoded_seq = self._kmer_one(stripped_seq)
        else:
            encoded_seq = self._kmer_n(stripped_seq)

        if seq[0] == self.START:
            encoded_seq = [self.encoding_index[self.START]] + encoded_seq

        non_seq_ending = "".join(s for s in seq if s in {self.END, self.PADDING})
        if non_seq_ending:
            encoded_seq = encoded_seq + [self.encoding_index[s] for s in non_seq_ending]

        return encoded_seq

    def integer_decode(self, int_seq: Sequence[int]) -> str:
        """Given a sequence of integers, convert it to a string."""

        return "".join(self.decode_token(s) for s in int_seq)
