# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""The base tokenizer from which other tokenizers are expected to inhereit.

The `SequenceTokenizerSettings` in particular holds common tokens and their associated ids for
common special tokens. For example `unk_token` for unknown tokens.
"""

from typing import Dict
from typing import List
from typing import Optional
from typing import Union

from pydantic import BaseSettings
from pydantic import Field
from pydantic import root_validator

DEFAULT_TOKEN_AND_IDS = {
    "pad": {"token": "|", "token_id": 0},
    "bos": {"token": ">", "token_id": 1},
    "eos": {"token": "<", "token_id": 2},
    "unk": {"token": "?", "token_id": 3},
    "mask": {"token": "#", "token_id": 4},
}


# pylint: disable=too-few-public-methods
class SequenceTokenizerSettings(BaseSettings):
    """The base tokenizer settings which inherits from pydantic.BaseSettings.

    To see the available fields, use pydantic's `schema_json`.

    ```python
    >>> print(SequenceTokenizerSettings.schema_json(indent=2))
    {
      "title": "SequenceTokenizerSettings"
      ...
    }
    """

    unk_token: Optional[str] = Field(None, env="GCGC_UNK_TOKEN")
    unk_token_id: Optional[int] = Field(None, env="GCGC_UNK_TOKEN_ID")

    bos_token: Optional[str] = Field(None, env="GCGC_BOS_TOKEN")
    bos_token_id: Optional[int] = Field(None, env="GCGC_BOS_TOKEN_ID")

    mask_token: Optional[str] = Field(None, env="GCGC_MASK_TOKEN")
    mask_token_id: Optional[int] = Field(None, env="GCGC_MASK_TOKEN_ID")

    eos_token: Optional[str] = Field(None, env="GCGC_EOS_TOKEN")
    eos_token_id: Optional[int] = Field(None, env="GCGC_EOS_TOKEN_ID")

    pad_token: Optional[str] = Field(None, env="GCGC_PAD_TOKEN")
    pad_token_id: Optional[int] = Field(None, env="GCGC_PAD_TOKEN_ID")

    max_length: Optional[int] = Field(None, env="GCGC_MAX_LENGTH")
    min_length: Optional[int] = Field(None, env="GCGC_MIN_LENGTH")

    conform_length: Optional[int] = Field(None, env="GCGC_CONFORM_LENGTH")

    def special_token_defaults(self) -> Dict[str, Dict[str, Union[int, str]]]:
        """Return a dictionary mapping the token to its default token and id."""
        return DEFAULT_TOKEN_AND_IDS

    # pylint: disable=no-self-argument, no-self-use
    @root_validator
    def validate_and_set_special_tokens(cls, values):
        """Check that the tokens and token id pairs have corresponding values.

        If they do not, this will set the default value.

        """
        passed_tokens = {
            k: v for k, v in values.items() if k.endswith("_token_id") and v is not None
        }

        for token_type, token_map in DEFAULT_TOKEN_AND_IDS.items():
            token_name = f"{token_type}_token"

            passed_token = values.get(token_name)
            passed_token_id = values.get(f"{token_type}_token_id")

            if (passed_token is None and passed_token_id is None) or (
                passed_token is not None and passed_token_id is not None
            ):
                continue

            if passed_token is None and passed_token_id is not None:
                values[token_name] = token_map["token"]

            if passed_token is not None and passed_token_id is None:
                if token_map["token_id"] in passed_tokens.values():
                    msg = (
                        f"Tried to set the default for {token_type}, however it "
                        f"conflicts with a manually supplied token_id in {passed_tokens}. "
                        "Supply an associated token value."
                    )
                    raise ValueError(msg)

                values[f"{token_name}_id"] = token_map["token_id"]

        return values

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
        """Init the sequence tokenizer with a set of settings."""
        self.settings = settings or SequenceTokenizerSettings()

    def get_special_tokens_mask(self, token_ids: List[int]) -> List[int]:
        """Given the input set of tokens, return a list that demarcates special character.

        ```python
        # Assuming 1 is bos_token, 2 is eos_token, and 0 is pad_token.
        >>> tokenizer.get_special_tokens_mask([1, 34, 21, 0, 0, 0, 2])
        [1, 0, 0, 1, 1, 1, 1]
        ```

        Args:
            token_ids: The list of integer tokens that may or may not be special tokens according
                to the toeknizer's settings.

        Returns:
            A list of 0s and 1s, where the value is one if the token is a special character.

        """
        return [int(token_id in self.settings.special_token_ids) for token_id in token_ids]

    def apply_length_constraints(self, tokens: List[str]) -> List[str]:
        """Apply the constraints from the settings to the passed tokens list."""
        if self.settings.conform_length:
            tokens = _pad_token_list(tokens, self.settings.conform_length, self.settings.pad_token)
            return tokens[: self.settings.conform_length]

        if self.settings.min_length:
            tokens = _pad_token_list(tokens, self.settings.min_length, self.settings.pad_token)

        if self.settings.max_length:
            tokens = tokens[: self.settings.max_length]

        return tokens
