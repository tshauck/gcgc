# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Modules for huggingface's transformers."""

from typing import List

from gcgc.tokenizer.kmer_tokenzier import KmerTokenizer

try:
    from tokenizers.pre_tokenizers import PreTokenizer

    from tokenizers import NormalizedString
    from tokenizers import PreTokenizedString
    from tokenizers import Tokenizer
    from tokenizers import models
    from tokenizers.decoders import ByteLevel

except ImportError as exp:
    # pylint: disable=invalid-name
    needed = "tokenizers"
    raise ImportError(f"Missing one or more libraries: {needed}. Please install: {exp}") from exp

# pylint: disable=too-few-public-methods


class KmerPreTokenizer:
    """Pretokenizes sequences based on kmers."""

    def __init__(self, kmer_length: int, kmer_stride: int, alphabet: str, unk_token: str = "?"):
        """Inits the KmerTokenizer.

        Args:
            kmer_length: How long of kmers to create.
            kmer_stride: The stride between two kmers. Should be equal to kmer_stride to generate
                non-overlapping kmers.
            alphabet: The particular alphabet
            unk_token: The unknown token to use for the pre-tokenization.

        """
        self.kmer_tokenzier = KmerTokenizer(
            kmer_length=kmer_length,
            kmer_stride=kmer_stride,
            alphabet=alphabet,
            pad_token=None,
            bos_token=None,
            eos_token=None,
            mask_token=None,
            unk_token=unk_token,
        )

    def pre_tokenize(self, pre_tok: PreTokenizedString):
        """Pretokenize the input string."""
        pre_tok.split(self._split)

    # pylint: disable=unused-argument
    def _split(self, i: int, normalized_string: NormalizedString) -> List[NormalizedString]:
        pre_tokenized_sequence = self.kmer_tokenzier.tokenize(str(normalized_string))
        return [NormalizedString(split) for split in pre_tokenized_sequence]


def build_hf_tokenizer(
    kmer_length: int, kmer_stride: int, alphabet: str, unk_token: str = "?"
) -> Tokenizer:
    """Build a full huggingface tokenizer from the inputs.

    Note:
        Same arguments taken as KmerPreTokenizer.

    """
    kmer_pre = KmerPreTokenizer(
        kmer_length=kmer_length, kmer_stride=kmer_stride, alphabet=alphabet, unk_token=unk_token
    )
    tokenizer = Tokenizer(
        models.WordLevel(vocab=kmer_pre.kmer_tokenzier.vocab.stoi, unk_token=unk_token)
    )
    tokenizer.pre_tokenizer = PreTokenizer.custom(kmer_pre)
    tokenizer.decoder = ByteLevel()

    return tokenizer
