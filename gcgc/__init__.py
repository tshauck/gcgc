# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Top-level GCGC module."""

import sys
import warnings as _warnings

from gcgc import tokenizer
from gcgc.tokenizer.kmer_tokenzier import KmerTokenizer
from gcgc.tokenizer.kmer_tokenzier import SequenceTokenizer

if sys.version_info[:2] >= (3, 8):
    # pylint: disable=no-name-in-module
    from importlib import metadata as importlib_metadata
else:
    import importlib_metadata

__version__ = importlib_metadata.version("gcgc")

__all__ = [
    "tokenizer",
    "KmerTokenizer",
    "SequenceTokenizer",
    "__version__",
]

_warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
