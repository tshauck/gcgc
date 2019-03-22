# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Utils for working with alphabets."""

from gcgc.alphabet import iupac
from gcgc.exceptions import GCGCAlphabetException


def biopython_alphabet_to_gcgc_alphabet(biopython_alphabet_instance, *args, **kwargs):
    """Convert the BioPython Alphabet to the associated GCGC Alphabet."""

    biopython_type = type(biopython_alphabet_instance)
    try:
        klass = iupac._ALPHABET_MAPPING[biopython_type]
        return klass(*args, **kwargs)
    except KeyError:
        msg = f"No instance of {biopython_type} among {list(iupac._ALPHABET_MAPPING)}."
        raise GCGCAlphabetException(msg)
