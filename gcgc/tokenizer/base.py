# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A Tokenizer that works on biological sequences."""

from typing import List, Optional

import pydantic


class SequenceTokenizerSpec(pydantic.BaseModel):
    """The specification for the tokenizer."""

    max_length: int = pydantic.Field(..., description="The maximum length of the sequence.")
    kmer_size: int = pydantic.Field(
        1, description="How big of kmers to tokenize the sequence into."
    )
    kmer_step_size: int = pydantic.Field(1, description="The size of the sliding window of kmer.")
    bos_token: Optional[str] = pydantic.Field(
        None, description="The token for the beginning of the sequence."
    )
    eos_token: Optional[str] = pydantic.Field(
        None, description="The token for the end of the sequence."
    )
    unk_token: Optional[str] = pydantic.Field(None, description="The token for an unknown token.")
    pad_token: Optional[str] = pydantic.Field(
        None, description="The token for the padding character."
    )
    cls_token: Optional[str] = pydantic.Field(None, description="The classification token.")
    mask_token: Optional[str] = pydantic.Field(
        None, description="The token for the mask character."
    )

    @property
    def passed_bos_token(self) -> bool:
        """Return True if this token is in use."""
        return self.bos_token is None

    @property
    def passed_eos_token(self) -> bool:
        """Return True if this token is in use."""
        return self.eos_token is None

    @property
    def passed_unk_token(self) -> bool:
        """Return True if this token is in use."""
        return self.unk_token is None

    @property
    def passed_pad_token(self) -> bool:
        """Return True if this token is in use."""
        return self.pad_token is None

    @property
    def passed_cls_token(self) -> bool:
        """Return True if this token is in use."""
        return self.cls_token is None

    @property
    def passed_mask_token(self) -> bool:
        """Return True if this token is in use."""
        return self.mask_token is None

    @property
    def max_encoded_length(self) -> int:
        """Returns how long the encoded sequence can be given the token's used."""
        return self.max_length - (int(self.passed_bos_token) + int(self.passed_eos_token))


class SequenceTokenizer:
    """A sequence tokenizer."""

    def __init__(self, tokenizer_spec: SequenceTokenizerSpec):
        """Init the SequenceTokenizer class.

        Args:
            tokenizer_spec: The spec for the tokenizer.

        """
        self.tokenizer_spec = tokenizer_spec

    @staticmethod
    def _kmer_one(seq: str) -> List[str]:
        return list(seq)

    def _kmer_n(self, seq: str) -> List[str]:
        seq_len = len(seq)
        iterations = seq_len - self.tokenizer_spec.kmer_size + 1
        kmers = []

        for i in range(0, iterations, self.tokenizer_spec.kmer_step_size):
            kmer = seq[i : i + self.tokenizer_spec.kmer_size]
            kmers.append(kmer)

        return kmers

    def __call__(self, seq: str) -> List[str]:
        """Integer encode the sequence.

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
            return self._kmer_one(seq)

        return self._kmer_n(seq)
