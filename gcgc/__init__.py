# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Top-level GCGC module."""

import os
import warnings as _warnings

from gcgc import tokenizer
from gcgc.tokenizer.kmer_tokenzier import KmerTokenizer
from gcgc.tokenizer.kmer_tokenzier import SequenceTokenizer

__all__ = ["tokenizer", "KmerTokenizer", "SequenceTokenizer"]

__version__ = "0.12.5"

_warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
