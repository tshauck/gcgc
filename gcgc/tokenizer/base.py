# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A Tokenizer that works on biological sequences."""

from pydantic import BaseSettings


class SequenceTokenizerSettings(BaseSettings):
    """The base tokenizer settings."""


class SequenceTokenizer:
    """The sequence tokenizer object."""
