# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A Tokenizer that works on biological sequences."""

from typing import List, Optional, Dict, cast
from dataclasses import field
import itertools as it

from pydantic import dataclasses


@dataclasses.dataclass
class Vocab:
    """A vocabulary object."""

    token_to_int: Dict[str, int]
    int_to_token: Dict[int, str]

    def __len__(self) -> int:
        """Return the length of the vocab."""
        return len(self.token_to_int)

    @classmethod
    def from_list(cls, tokens: List[str]) -> "Vocab":
        """Create a vocabulary from a list of tokens."""
        token_to_int = {}
        int_to_token = {}

        for i, token in enumerate(tokens):
            int_to_token[i] = token
            token_to_int[token] = i

        return cls(token_to_int, int_to_token)


@dataclasses.dataclass
class SequenceTokenizerSpec:
    """The specification for the tokenizer."""

    max_length: int
    alphabet: str

    vocabulary: Vocab = field(init=False)

    kmer_size: int = 1
    kmer_step_size: int = 1

    bos_token: Optional[str] = None
    eos_token: Optional[str] = None
    unk_token: Optional[str] = None
    pad_token: Optional[str] = None
    mask_token: Optional[str] = None

    def __post_init__(self):
        """Post inits the tokenizer spec."""
        self.vocabulary = Vocab.from_list(self.special_tokens + self.possible_kmers)

    @property
    def possible_kmers(self) -> List[str]:
        """Return the set of possible kmers given the alphabet and kmer size."""
        return ["".join(kmer) for kmer in it.product(self.alphabet, repeat=self.kmer_size)]

    @property
    def special_tokens(self) -> List[str]:
        """Return the list of special tokens."""
        return [
            s
            for s in [
                self.bos_token,
                self.eos_token,
                self.unk_token,
                self.pad_token,
                self.mask_token,
            ]
            if s
        ]

    @property
    def passed_bos_token(self) -> bool:
        """Return True if this token is in use."""
        return self.bos_token is not None

    @property
    def passed_eos_token(self) -> bool:
        """Return True if this token is in use."""
        return self.eos_token is not None

    @property
    def passed_unk_token(self) -> bool:
        """Return True if this token is in use."""
        return self.unk_token is not None

    @property
    def passed_pad_token(self) -> bool:
        """Return True if this token is in use."""
        return self.pad_token is not None

    @property
    def passed_mask_token(self) -> bool:
        """Return True if this token is in use."""
        return self.mask_token is not None

    @property
    def max_tokenized_length(self) -> int:
        """Return how long the tokenized list can be given the tokens used."""
        return self.max_length - (int(self.passed_bos_token) + int(self.passed_eos_token))

    @property
    def bos_token_int(self) -> int:
        """Return the integer encoding of the passed token."""
        if not self.passed_bos_token:
            raise ValueError(f"bos_token is false-y ({self.bos_token}), cannot get int.")
        return self.vocabulary.token_to_int[cast(str, self.bos_token)]

    @property
    def eos_token_int(self) -> int:
        """Return the integer encoding of the passed token."""
        if not self.passed_eos_token:
            raise ValueError(f"eos_token is false-y ({self.eos_token}), cannot get int.")
        return self.vocabulary.token_to_int[cast(str, self.eos_token)]

    @property
    def unk_token_int(self) -> int:
        """Return the integer encoding of the passed token."""
        if not self.passed_unk_token:
            raise ValueError(f"unk_token is false-y ({self.unk_token}), cannot get int.")
        return self.vocabulary.token_to_int[cast(str, self.unk_token)]

    @property
    def pad_token_int(self) -> int:
        """Return the integer encoding of the passed token."""
        if not self.passed_pad_token:
            raise ValueError(f"pad_token is false-y ({self.pad_token}), cannot get int.")
        return self.vocabulary.token_to_int[cast(str, self.pad_token)]

    @property
    def mask_token_int(self) -> int:
        """Return the integer encoding of the passed token."""
        if not self.passed_mask_token:
            raise ValueError(f"mask_token is false-y ({self.mask_token}), cannot get int.")
        return self.vocabulary.token_to_int[cast(str, self.mask_token)]


class SequenceTokenizer:
    """A sequence tokenizer."""

    def __init__(self, tokenizer_spec: SequenceTokenizerSpec):
        """Init the SequenceTokenizer class.

        Args:
            tokenizer_spec: The spec for the tokenizer.

        """
        self.tokenizer_spec = tokenizer_spec

    def _kmer_n(self, seq: str) -> List[str]:
        seq_len = len(seq)
        iterations = seq_len - self.tokenizer_spec.kmer_size + 1
        kmers = []

        for i in range(0, iterations, self.tokenizer_spec.kmer_step_size):
            kmer = seq[i : i + self.tokenizer_spec.kmer_size]
            kmers.append(kmer)

        return kmers

    def __call__(self, seq: str) -> List[str]:
        """Tokenize the seqeunce, see .tokenize."""
        return self.tokenize(seq)

    def encode(self, seq: str) -> List[int]:
        """Encode the underlying sequence."""
        return [self.tokenizer_spec.vocabulary.token_to_int[s] for s in self.tokenize(seq)]

    def tokenize(self, seq: str) -> List[str]:
        """Tokenize the sequence.

        Args:
            seq: The sequence to encode.

        Returns:
            The list of strs that are the tokens.

        """
        seq_len = len(seq)

        if seq_len < self.tokenizer_spec.kmer_size:
            raise ValueError(
                f"seq length {seq_len} cannot be less than the kmer "
                "size {self.tokenizer_spec.kmer_size}"
            )

        if self.tokenizer_spec.kmer_size == 1:
            kmer_list = list(seq)
        else:
            kmer_list = self._kmer_n(seq)

        sequence_kmers = kmer_list[: self.tokenizer_spec.max_tokenized_length]

        if self.tokenizer_spec.passed_bos_token:
            sequence_kmers.insert(0, cast(str, self.tokenizer_spec.bos_token))

        if self.tokenizer_spec.passed_eos_token:
            sequence_kmers.append(cast(str, self.tokenizer_spec.eos_token))

        return sequence_kmers
