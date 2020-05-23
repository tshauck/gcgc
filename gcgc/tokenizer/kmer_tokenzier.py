# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""This module holds the kmer tokenizer and its settings.

The `KmerTokenizer` tokenizers incoming sequences into kmers of configurable size against a
configurable alphabet.

For example, given incoming nucleotide sequences, using settings of kmer length 3 and stride 3 they
encoded sequences will be codons (in the loose sense) with a vocabulary of size 64.

```python

>>> tokenizer.encode('AAATTTCCC')
[0, 42, 63]

>>> len(tokenizer.vocab)
64
```
"""

import itertools as it
from typing import List

from pydantic import Field
from pydantic import root_validator
from pydantic import validator

from gcgc import alphabets
from gcgc.tokenizer.base import SequenceTokenizer
from gcgc.vocab import Vocab


def _create_kmer_vocab_from_token(vocab, alphabet: str, kmer_length: int) -> Vocab:
    """Create the vocab object from a list of tokens."""

    token_set = ["".join(kmer) for kmer in it.product(list(alphabet), repeat=kmer_length)]

    for token in token_set:
        vocab.add_item(token)

    return vocab


class KmerTokenizer(SequenceTokenizer):
    """The Kmer Tokenizer that encodes sequences into chunked kmers."""

    alphabet: str = Field("ATCG", env="GCGC_ALPHABET")
    kmer_length: int = Field(1, env="GCGC_KMER_LENGTH")
    kmer_stride: int = Field(1, env="GCGC_KMER_STRIDE")

    @root_validator
    def augment_vocab_with_kmer(cls, values):  # pylint: disable=no-self-argument,no-self-use
        """Update the vocab of the SequenceTokenizer with a kmer alphabet."""
        vocab = values["vocab"]
        alphabet = values["alphabet"]
        kmer_length = values["kmer_length"]

        values["vocab"] = _create_kmer_vocab_from_token(vocab, alphabet, kmer_length)
        return values

    # pylint: disable=no-self-use, no-self-argument
    @validator("alphabet")
    def resolve_alphabet(cls, alphabet: str) -> str:
        """Resolve the alphabet if it's a named alphabet.

        Args:
            alphabet: The raw alphabet, either the sequence literal or a name of a preset alphabet.

        Returns:
            The new alphabet.

        """
        return alphabets.resolve_alphabet(alphabet)

    def _kmer_n(self, seq: str) -> List[str]:
        seq_len = len(seq)
        iterations = seq_len - self.kmer_length + 1
        kmers = []

        for i in range(0, iterations, self.kmer_stride):
            kmer = seq[i : i + self.kmer_length]
            kmers.append(kmer)

        return kmers

    def encode_batches(self, seqs: List[str], add_unknown: bool = True) -> List[List[int]]:
        """Encode batches rather than a single example.

        Args:
            seqs: A list of sequences.
            add_unknown: Passed to encode.

        Returns:
            A list of encoded sequences.

        """
        return [self.encode(seq, add_unknown=add_unknown) for seq in seqs]

    def encode(self, seq: str, add_unknown: bool = False) -> List[int]:
        """Encode the underlying sequence into a list of tokens ids.

        Args:
            seq: The incoming sequence to encode.
            add_unknown: If True, add the unknown token rather than throwing an out of vocabulary
                error.

        Returns:
            A list of the encoded tokens.

        """
        encoded = []

        for letter in self.encode_as_tokens(seq):
            try:
                encoded.append(self.vocab[letter])
            except KeyError:
                if add_unknown:
                    encoded.append(self.unk_token_id)
                else:
                    raise

        return encoded

    def encode_as_tokens(self, seq: str) -> List[str]:
        """Tokenize the sequence into a list of tokens.

        Args:
            seq: The sequence to encode.

        Returns:
            The list of strings that are the tokens.

        """
        seq_len = len(seq)

        if seq_len < self.kmer_length:
            raise ValueError(
                f"seq length {seq_len} cannot be less than the kmer size {self.kmer_length}"
            )

        if self.kmer_length == 1:
            kmer_list = list(seq)
        else:
            kmer_list = self._kmer_n(seq)

        if self.bos_token:
            kmer_list = [self.bos_token] + kmer_list

        if self.eos_token:
            kmer_list = kmer_list + [self.eos_token]

        return super().apply_length_constraints(kmer_list)
