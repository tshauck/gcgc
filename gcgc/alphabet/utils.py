# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Utils for working with alphabets."""

from gcgc import alphabet
from gcgc.exceptions import GCGCAlphabetException


def biopython_alphabet_to_gcgc_alphabet(biopython_alphabet_instance):
    """Convert the BioPython Alphabet to the associated GCGC Alphabet."""

    # Might also try just creating the type on the fly.

    for gcgc_alphabet in alphabet.__all__:

        klass = getattr(alphabet, gcgc_alphabet)
        instance = klass()

        if isinstance(instance, type(biopython_alphabet_instance)):
            return instance

    msg = f"No instance of {type(biopython_alphabet_instance)} among {alphabet.__all__}."
    raise GCGCAlphabetException(msg)
