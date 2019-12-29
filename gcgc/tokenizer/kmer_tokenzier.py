# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Module for KMER tokenization."""

from typing import List, Optional
import itertools as it

from pydantic import Field, validator

from gcgc.tokenizer.base import SequenceTokenizer, SequenceTokenizerSettings
from gcgc import alphabets


class KmerTokenizerSettings(SequenceTokenizerSettings):
    """The specification for the tokenizer."""

    alphabet: str = Field("ATCG", env="GCGC_ALPHABET")
    kmer_length: int = Field(1, env="GCGC_KMER_LENGTH")
    kmer_stride: int = Field(1, env="GCGC_KMER_STRIDE")

    @validator("alphabet")
    def resolve_alphabet(cls, alphabet):  # pylint: disable=no-self-use, no-self-argument
        """Resolve the alphabet if it's a named alphabet."""
        return alphabets.resolve_alphabet(alphabet)

    @property
    def special_tokens(self) -> List[str]:
        """Return the list of special tokens that are not None."""
        return [
            x
            for x in [
                self.bos_token,
                self.eos_token,
                self.unk_token,
                self.pad_token,
                self.mask_token,
            ]
            if x is not None
        ]


def _create_kmer_vocab_from_token(
    alphabet, kmer_length, token_characters: Optional[List[str]] = None
):
    """Create vocabulary object from a list of tokens."""
    if token_characters is None:
        token_characters = []

    token_to_int = {}

    for i, token in enumerate(token_characters):
        token_to_int[token] = i

    token_list = ["".join(kmer) for kmer in it.product(list(alphabet), repeat=kmer_length)]
    for i, token in enumerate(token_list, start=len(token_characters)):
        token_to_int[token] = i

    # pylint: disable=too-many-function-args
    return token_to_int


class KmerTokenizer(SequenceTokenizer):
    """A sequence tokenizer."""

    def __init__(self, settings: Optional[KmerTokenizerSettings] = None):
        """Init the SequenceTokenizer class.

        Args:
            settings: The settings for the tokenizer.
            vocabulary: The vocabulary for the tokenizer.

        """
        self.settings = settings or KmerTokenizerSettings()
        super().__init__(settings)

        self.vocab = _create_kmer_vocab_from_token(
            self.settings.alphabet, self.settings.kmer_length, self.settings.special_tokens,
        )

    def _kmer_n(self, seq: str) -> List[str]:
        seq_len = len(seq)
        iterations = seq_len - self.settings.kmer_length + 1
        kmers = []

        for i in range(0, iterations, self.settings.kmer_stride):
            kmer = seq[i : i + self.settings.kmer_length]
            kmers.append(kmer)

        return kmers

    def encode(self, seq: str) -> List[int]:
        """Encode the underlying sequence into a list of tokens."""
        return [
            self.vocab.get(s, self.vocab.get(self.settings.unk_token))
            for s in self.encode_as_tokens(seq)
        ]

    def encode_as_tokens(self, seq: str) -> List[str]:
        """Tokenize the sequence into a list of token tokens.

        Args:
            seq: The sequence to encode.

        Returns:
            The list of strs that are the tokens.

        """
        seq_len = len(seq)

        if seq_len < self.settings.kmer_length:
            raise ValueError(
                f"seq length {seq_len} cannot be less than the kmer "
                "size {self.settings.kmer_length}"
            )

        if self.settings.kmer_length == 1:
            kmer_list = list(seq)
        else:
            kmer_list = self._kmer_n(seq)

        if self.settings.bos_token:
            kmer_list = [self.settings.bos_token] + kmer_list

        if self.settings.eos_token:
            kmer_list = kmer_list + [self.settings.eos_token]

        return super().apply_length_constraints(kmer_list)
