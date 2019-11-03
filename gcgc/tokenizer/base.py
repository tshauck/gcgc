# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A Tokenizer that works on biological sequences."""

from typing import List


class SequenceTokenizer:
    """A sequence tokenizer."""

    def __init__(self, kmer_size: int = 1, kmer_step_size: int = 1):
        """Init the SequenceTokenizer class.

        Args:
            kmer_size: How big of kmers to tokenize the sequence into.

        """
        self.kmer_size = kmer_size
        self.kmer_step_size = kmer_step_size

    @staticmethod
    def _kmer_one(seq: str) -> List[str]:
        return list(seq)

    def _kmer_n(self, seq: str) -> List[str]:
        seq_len = len(seq)
        iterations = seq_len - self.kmer_size + 1
        kmers = []

        for i in range(0, iterations, self.kmer_step_size):
            kmer = seq[i : i + self.kmer_size]
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

        if seq_len < self.kmer_size:
            raise ValueError(
                f"seq length {seq_len} cannot be less than the kmer size {self.kmer_size}"
            )

        if self.kmer_size == 1:
            return self._kmer_one(seq)

        return self._kmer_n(seq)
