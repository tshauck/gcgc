# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""The base tokenizer."""

from typing import List
from typing import Optional

from pydantic import BaseSettings
from pydantic import Field
from pydantic import root_validator


# pylint: disable=too-few-public-methods
class SequenceTokenizerSettings(BaseSettings):
    """The base tokenizer settings."""

    bos_token: Optional[str] = Field(None, env="GCGC_BOS_TOKEN")
    eos_token: Optional[str] = Field(None, env="GCGC_EOS_TOKEN")
    unk_token: Optional[str] = Field("?", env="GCGC_UNK_TOKEN")
    mask_token: Optional[str] = Field(None, env="GCGC_MASK_TOKEN")
    pad_token: Optional[str] = Field(None, env="GCGC_PAD_TOKEN")

    max_length: Optional[int] = Field(None, env="GCGC_MAX_LENGTH")
    min_length: Optional[int] = Field(None, env="GCGC_MIN_LENGTH")
    conform_length: Optional[int] = Field(None, env="GCGC_CONFORM_LENGTH")

    @root_validator
    def validate_lengths(cls, values):  # pylint: disable=no-self-argument, no-self-use
        """Check the length arguments are valid."""
        max_length = values.get("max_length")
        min_length = values.get("min_length")
        conform_length = values.get("conform_length")
        pad_token = values.get("pad_token")

        if conform_length is not None and (max_length is not None or min_length is not None):
            raise ValueError(f"If conform length is not None, max and min length can't be set.")

        if max_length is not None and min_length is not None and min_length > max_length:
            raise ValueError(
                f"Min length {min_length} cannot be more than max length {max_length}."
            )

        if (conform_length is not None or min_length is not None) and pad_token is None:
            raise ValueError("Cannot set conform_length or min_length without a pad_token.")

        return values


def _pad_token_list(tokens: List[str], pad_to: int, pad_char: str):
    """Pad the input tokens list to tokens."""
    token_len = len(tokens)

    if token_len > pad_to:
        return tokens

    new_tokens = pad_to - token_len
    return tokens + [pad_char] * new_tokens


class SequenceTokenizer:
    """The sequence tokenizer object."""

    def __init__(self, settings: Optional[SequenceTokenizerSettings] = None):
        """Inits the sequence tokenizer with a set of settings."""
        self.settings = settings or SequenceTokenizerSettings()

    def apply_length_constraints(self, tokens: List[str]):
        """Apply the constraints of the settings to the passed tokens list."""
        if self.settings.conform_length:
            tokens = _pad_token_list(tokens, self.settings.conform_length, self.settings.pad_token)
            return tokens[: self.settings.conform_length]

        if self.settings.min_length:
            tokens = _pad_token_list(tokens, self.settings.min_length, self.settings.pad_token)

        if self.settings.max_length:
            tokens = tokens[: self.settings.max_length]

        return tokens
