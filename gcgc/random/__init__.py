# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Module to help generate random sequences."""

import numpy as np
from gcgc.alphabet.base import EncodingAlphabet
from gcgc.encoded_seq import EncodedSeq


def random_sequence_from_alphabet(alphabet: EncodingAlphabet, n: int) -> EncodedSeq:
    """Randomly sample from the alphabet to length n."""

    random_letters = np.random.choice(list(alphabet.letters_and_tokens), n)
    return EncodedSeq("".join(random_letters), alphabet)
