# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""GCGC Exceptions."""


class GCGCAlphabetException(Exception):
    """Raised when there's an error with GCGC's alphabet."""


class GCGCAlphabetLetterEncodingException(GCGCAlphabetException):
    """Thrown when GCGC's alphabet tries to encode a letter it wasn't aware of."""


class GCGCAlphabetLetterDecodingException(GCGCAlphabetException):
    """Thrown when GCGC's alphabet tries to decode a letter it wasn't aware of."""


class EncodedSeqLengthParserException(Exception):
    """Raised when there's an issue with the EncodedSeqLengthParser's parsing."""
