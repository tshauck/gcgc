# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Module for KMER tokenization."""

from typing import List, Optional
import itertools as it
from pathlib import Path

from pydantic import Field

from gcgc.tokenizer.base import SequenceTokenizer, SequenceTokenizerSettings

try:
    import sentencepiece as spm

    has_spm = True
except ImportError:
    has_spm = False


class BioSequencePieceSettings(SequenceTokenizerSettings):
    """The settings for the sentence piece model."""

    model_path: Path = Field(..., env="GCGC_SP_MODEL_PATH")
    vocab_path: Path = Field(..., env="GCGC_SP_VOCAB_PATH")

    max_length: Optional[int] = Field(None, env="GCGC_KMER_MAX_LENGTH")

    bos_token: Optional[str] = Field(None, env="GCGC_BOS_TOKEN")
    eos_token: Optional[str] = Field(None, env="GCGC_EOS_TOKEN")
    unk_token: Optional[str] = Field(None, env="GCGC_UNK_TOKEN")
    pad_token: Optional[str] = Field(None, env="GCGC_PAD_TOKEN")
    mask_token: Optional[str] = Field(None, env="GCGC_MASK_TOKEN")

    @property
    def special_tokens(self) -> List[str]:
        """Return the special tokens that are not None."""
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


class BioSequencePiece(SequenceTokenizer):
    """A sentence piece for model on biological sequences."""

    def __init__(self, tokenizer_settings: Optional[BioSequencePieceSettings] = None):
        """Init the BioSequencePiece class.

        Args:
            tokenizer_settings: The settings for the tokenizer.
            vocabulary: The vocabulary for the tokenizer.

        """
        if not has_spm:
            raise RuntimeError(f"Trying to use sentencepiece but the python library is missing!")

        self.tokenizer_settings = tokenizer_settings or BioSequencePieceSettings()
        self.sp_processor = spm.SentencePieceProcessor()
        self.sp_processor.load(self.tokenizer_settings.model_path)

    def encode_as_ids(self, seq: str) -> List[int]:
        """Encode the underlying sequence into a list of tokens."""
        return self.sp_processor.EncodeAsIds(seq)

    def encode_as_tokens(self, seq: str) -> List[str]:
        """Tokenize the sequence into a list of token tokens.

        Args:
            seq: The sequence to encode.

        Returns:
            The list of strs that are the tokens.

        """
        return self.sp_processor.EncodeAsPieces(seq)
