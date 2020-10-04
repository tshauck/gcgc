# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""The base tokenizer from which other tokenizers are expected to inhereit.

The `SequenceTokenizerSettings` in particular holds common tokens and their associated ids for
common special tokens. For example `unk_token` for unknown tokens.
"""

from typing import List
from typing import Optional

from pydantic import BaseSettings
from pydantic import Field
from pydantic import root_validator

from gcgc.vocab import Vocab


def _pad_token_list(tokens: List[str], pad_to: int, pad_char: str, pad_at_end: bool):
    """Pad the input tokens list to tokens."""
    token_len = len(tokens)

    if token_len > pad_to:
        return tokens

    new_tokens = pad_to - token_len

    padding = [pad_char] * new_tokens
    if pad_at_end:
        return tokens + padding
    return padding + tokens


class SequenceTokenizer(BaseSettings):
    """The sequence tokenizer object."""

    vocab: Vocab = Vocab()

    pad_token: Optional[str] = Field("|", env="GCGC_PAD_TOKEN")
    bos_token: Optional[str] = Field(">", env="GCGC_BOS_TOKEN")
    eos_token: Optional[str] = Field("<", env="GCGC_EOS_TOKEN")
    mask_token: Optional[str] = Field("#", env="GCGC_MASK_TOKEN")
    unk_token: Optional[str] = Field("?", env="GCGC_UNK_TOKEN")

    pad_token_id: Optional[int] = Field(None, env="GCGC_PAD_TOKEN_ID")
    bos_token_id: Optional[int] = Field(None, env="GCGC_BOS_TOKEN_ID")
    eos_token_id: Optional[int] = Field(None, env="GCGC_EOS_TOKEN_ID")
    mask_token_id: Optional[int] = Field(None, env="GCGC_MASK_TOKEN_ID")
    unk_token_id: Optional[int] = Field(None, env="GCGC_UNK_TOKEN_ID")

    pad_at_end: bool = Field(True, env="GCGC_PAD_AT_END")

    max_length: Optional[int] = Field(None, env="GCGC_MAX_LENGTH")
    min_length: Optional[int] = Field(None, env="GCGC_MIN_LENGTH")
    conform_length: Optional[int] = Field(None, env="GCGC_CONFORM_LENGTH")

    @classmethod
    def bare_tokenizer(cls, **values):
        """Init a tokenizer like normal, but default all tokens to None."""
        pad_token = values.pop("pad_token", None)
        bos_token = values.pop("bos_token", None)
        eos_token = values.pop("eos_token", None)
        mask_token = values.pop("mask_token", None)
        unk_token = values.pop("unk_token", None)

        return cls(
            pad_token=pad_token,
            bos_token=bos_token,
            eos_token=eos_token,
            mask_token=mask_token,
            unk_token=unk_token,
            **values,
        )

    @root_validator
    def token_ids(cls, values):  # pylint: disable=no-self-use,no-self-argument
        """Update the pad token id if appropriate."""
        vocab = values["vocab"]

        for key, value in values.items():
            if key.endswith("_token"):
                token_id_field = f"{key}_id"
                token_id = values[token_id_field]

                if token_id is None and value is not None:
                    vocab.add_item(value)
                    values[token_id_field] = vocab[value]
                elif token_id is not None and value is not None:
                    # pylint: disable=fixme
                    # TODO: check that token_id isn't already set (maybe its in the inverse?)
                    vocab[value] = token_id

        return values

    @root_validator
    def validate_lengths(cls, values):  # pylint: disable=no-self-argument, no-self-use
        """Check the length arguments are valid."""
        max_length = values.get("max_length")
        min_length = values.get("min_length")
        conform_length = values.get("conform_length")
        pad_token = values.get("pad_token")

        if conform_length is not None and (max_length is not None or min_length is not None):
            raise ValueError("If conform length is not None, max and min length can't be set.")

        if max_length is not None and min_length is not None and min_length > max_length:
            raise ValueError(
                f"Min length {min_length} cannot be more than max length {max_length}."
            )

        if (conform_length is not None or min_length is not None) and pad_token is None:
            raise ValueError("Cannot set conform_length or min_length without a pad_token.")

        return values

    @property
    def special_token_ids(self) -> List[int]:
        """Return the list of token ids corresponding to special tokens."""
        return [
            x
            for x in [
                self.bos_token_id,
                self.eos_token_id,
                self.unk_token_id,
                self.pad_token_id,
                self.mask_token_id,
            ]
            if x is not None
        ]

    @property
    def special_tokens(self) -> List[str]:
        """Return the set of special tokens that are not None."""
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

    def get_special_tokens_mask(self, token_ids: List[int]) -> List[int]:
        """Given the input set of tokens, return a list that demarcates special character.

            >>> from gcgc import KmerTokenizer
            >>> tokenizer = KmerTokenizer()

            # Assuming 1 is bos_token, 2 is eos_token, and 0 is pad_token.
            >>> tokenizer.get_special_tokens_mask([1, 34, 21, 0, 0, 0, 2])
            [1, 0, 0, 1, 1, 1, 1]

        Args:
            token_ids: The list of integer tokens that may or may not be special tokens according
                to the toeknizer's settings.

        Returns:
            A list of 0s and 1s, where the value is one if the token is a special character.

        """
        return [int(token_id in self.special_token_ids) for token_id in token_ids]

    def apply_length_constraints(self, tokens: List[str]) -> List[str]:
        """Apply the constraints from the settings to the passed tokens list."""
        if self.conform_length:
            tokens = _pad_token_list(
                tokens, self.conform_length, str(self.pad_token), self.pad_at_end,
            )
            return tokens[: self.conform_length]

        if self.min_length:
            tokens = _pad_token_list(tokens, self.min_length, str(self.pad_token), self.pad_at_end)

        if self.max_length:
            tokens = tokens[: self.max_length]

        return tokens
