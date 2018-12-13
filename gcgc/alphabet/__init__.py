# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""GCGC Alphabets."""

from gcgc.alphabet.iupac import ExtendedIUPACProteinEncoding
from gcgc.alphabet.iupac import ExtendedIUPACDNAEncoding
from gcgc.alphabet.iupac import IUPACUnambiguousDNAEncoding
from gcgc.alphabet.iupac import IUPACUnambiguousRNAEncoding
from gcgc.alphabet.iupac import IUPACAmbiguousDNAEncoding
from gcgc.alphabet.iupac import IUPACAmbiguousRNAEncoding
from gcgc.alphabet.iupac import IUPACProteinEncoding

__all__ = [
    "ExtendedIUPACProteinEncoding",
    "ExtendedIUPACDNAEncoding",
    "IUPACUnambiguousDNAEncoding",
    "IUPACUnambiguousRNAEncoding",
    "IUPACAmbiguousDNAEncoding",
    "IUPACAmbiguousRNAEncoding",
    "IUPACProteinEncoding",
]
