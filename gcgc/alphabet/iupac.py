# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from Bio.Alphabet import IUPAC

from gcgc.alphabet.base import EncodingAlphabet


class ExtendedIUPACProteinEncoding(EncodingAlphabet, IUPAC.ExtendedIUPACProtein):
    """
    Implements an encoding alphabet using the IUPAC extended protein letters.
    """


class ExtendedIUPACDNAEncoding(EncodingAlphabet, IUPAC.ExtendedIUPACDNA):
    """
    Implements an encoding alphabet using the IUPAC extended dna letters.
    """


class IUPACDnaAlphabetEncoding(EncodingAlphabet, IUPAC.IUPACUnambiguousDNA):
    """
    Implements an encoding alphabet using the IUPAC unambiguous dna letters.
    """
