# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from Bio.Data import IUPACData

from gcgc.alphabet.base import EncodingAlphabet


class ExtendedProteinAlphabet(EncodingAlphabet):
    @property
    def letters(self) -> str:
        return IUPACData.extended_protein_letters


class ExtendedDnaAlphabet(EncodingAlphabet):
    @property
    def letters(self) -> str:
        return IUPACData.extended_dna_letters


class UnambiguousDnaAlphabet(EncodingAlphabet):
    @property
    def letters(self) -> str:
        return IUPACData.unambiguous_dna_letters
