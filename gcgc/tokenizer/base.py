# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A Tokenizer that works on biological sequences."""

from typing import Optional

from pydantic import BaseSettings, Field


class SequenceTokenizerSettings(BaseSettings):
    """The base tokenizer settings."""

    bos_token: Optional[str] = Field(None, env="GCGC_BOS_TOKEN")
    eos_token: Optional[str] = Field(None, env="GCGC_EOS_TOKEN")
    unk_token: Optional[str] = Field("?", env="GCGC_UNK_TOKEN")
    mask_token: Optional[str] = Field(None, env="GCGC_MASK_TOKEN")
    pad_token: Optional[str] = Field(None, env="GCGC_PAD_TOKEN")


class SequenceTokenizer:
    """The sequence tokenizer object."""
