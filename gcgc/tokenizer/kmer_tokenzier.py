# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Module for KMER tokenization."""

import itertools as it
from typing import Dict
from typing import List
from typing import Optional

from pydantic import Field
from pydantic import validator

from gcgc import alphabets
from gcgc.tokenizer.base import SequenceTokenizer
from gcgc.tokenizer.base import SequenceTokenizerSettings


class KmerTokenizerSettings(SequenceTokenizerSettings):
    """The specification for the tokenizer."""

    alphabet: str = Field("ATCG", env="GCGC_ALPHABET")
    kmer_length: int = Field(1, env="GCGC_KMER_LENGTH")
    kmer_stride: int = Field(1, env="GCGC_KMER_STRIDE")

    @validator("alphabet")
    def resolve_alphabet(cls, alphabet):  # pylint: disable=no-self-use, no-self-argument
        """Resolve the alphabet if it's a named alphabet."""
        return alphabets.resolve_alphabet(alphabet)


def _create_kmer_vocab_from_token(settings: KmerTokenizerSettings) -> Dict[str, int]:
    """Create vocabulary object from a list of tokens."""
    token_to_int = {}

    for token_id, token in zip(settings.special_token_ids, settings.special_tokens):
        token_to_int[token] = token_id

    token_set = [
        "".join(kmer) for kmer in it.product(list(settings.alphabet), repeat=settings.kmer_length)
    ]

    starting_value = 0
    for token in token_set:
        while starting_value in settings.special_token_ids:
            starting_value += 1

        token_to_int[token] = starting_value
        starting_value += 1

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

        self.vocab = _create_kmer_vocab_from_token(self.settings)

    def _kmer_n(self, seq: str) -> List[str]:
        seq_len = len(seq)
        iterations = seq_len - self.settings.kmer_length + 1
        kmers = []

        for i in range(0, iterations, self.settings.kmer_stride):
            kmer = seq[i : i + self.settings.kmer_length]
            kmers.append(kmer)

        return kmers

    def encode(self, seq: str, add_unknown: bool = False) -> List[int]:
        """Encode the underlying sequence into a list of tokens."""
        encoded = []

        for letter in self.encode_as_tokens(seq):
            try:
                if self.vocab[letter] == 2 and letter == ">":
                    __import__("ipdb").set_trace()
                encoded.append(self.vocab[letter])
            except KeyError:
                if add_unknown:
                    encoded.append(self.settings.unk_token_id)
                else:
                    raise

        return encoded

        # return [
        # self.vocab.get(s, self.vocab.get(self.settings.unk_token))
        # for s in self.encode_as_tokens(seq)
        # ]

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
