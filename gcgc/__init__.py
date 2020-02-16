# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Top-level GCGC module."""

import warnings as _warnings

from gcgc import tokenizer
from gcgc.cli import cli
from gcgc.tokenizer import KmerTokenizer
from gcgc.tokenizer import KmerTokenizerSettings
from gcgc.tokenizer import SequenceTokenizer
from gcgc.tokenizer import SequenceTokenizerSettings

__version__ = "0.12.1-dev.3"

_warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
