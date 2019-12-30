# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Stores the alphabets, these are the same as BioPython."""

_ALPHABETS = {
    "protein": "ACDEFGHIKLMNPQRSTVWY",
    "extended_protein": "ACDEFGHIKLMNPQRSTVWYBXZJUO",
    "ambiguous_dna": "GATCRYWSMKHBVDN",
    "unambiguous_dna": "GATC",
    "extended_dna": "GATCBDSW",
    "ambiguous_rna": "GAUCRYWSMKHBVDN",
    "unambiguous_rna": "GAUC",
}


def resolve_alphabet(alphabet: str) -> str:
    """Try to get the alphabet from the known list, otherwise, return what was passed.

    Args:
        alphabet: The alphabet to resolve.

    Returns:
        The resolved alphabet.

    """
    try:
        return _ALPHABETS[alphabet]
    except KeyError:
        return alphabet
