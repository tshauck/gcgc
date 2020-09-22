# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Stores a common set of alphabets."""

ALPHABETS = {
    "protein": "ACDEFGHIKLMNPQRSTVWY",
    "extended_protein": "ACDEFGHIKLMNPQRSTVWYBXZJUO",
    "ambiguous_dna": "GATCRYWSMKHBVDN",
    "unambiguous_dna": "GATC",
    "extended_dna": "GATCBDSW",
    "ambiguous_rna": "GAUCRYWSMKHBVDN",
    "unambiguous_rna": "GAUC",
}


def resolve_alphabet(alphabet: str, require: bool = False) -> str:
    """Try to get the alphabet from the known list, otherwise, return what was passed.

    Args:
        alphabet: The alphabet to resolve.
        require: If true, require the passed alphabet be found in the options, else error.
            Defaults to false. This is useful if you want to use a named alphabet, but want to make
            sure it's valid.

    Returns:
        The resolved alphabet.

    """
    try:
        return ALPHABETS[alphabet]
    except KeyError:
        if require:
            raise
        return alphabet
