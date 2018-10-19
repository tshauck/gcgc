# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from typing import Sequence

from Bio.SeqRecord import SeqRecord
import numpy as np


class EncodedSeqRecord(object):
    """
    An inteface that holds information on encoded types.
    """

    def __init__(
        self,
        alphabet,
        seq_record: SeqRecord,
        encapsulate: bool = True,
        padding_to: int = 50,
    ) -> None:

        self.alphabet = alphabet
        self.seq_record = seq_record
        self.encapsulate = encapsulate
        self.padding_to = padding_to

    @property
    def _encapsulated_sequence(self) -> str:
        """
        Wrap a sequence with the start and end tokens.

        Returns:
            The sequence wrapped in START and END from the alphabet.
        """
        return f"{self.alphabet.START}{self.seq_record.seq}{self.alphabet.END}"

    @property
    def _seq(self) -> str:
        """
        Return the sequence that is encapsulated in the tokens, if the associated flag is set.
        """

        if self.encapsulate:
            return self._encapsulated_sequence
        else:
            return self.seq_record.seq

    @property
    def padded(self) -> str:
        """
        Pad a sequence up to self.padding_to if it's shorter, otherwise return the sequence.

        Returns:
            The sequence padding up to padding_to characters.
        """

        working_seq = self._seq
        seq_len = len(working_seq)

        if seq_len < self.padding_to:
            extra_chars = self.padding_to - seq_len
            encoding_seq = working_seq + (self.alphabet.PADDING * extra_chars)
        else:
            encoding_seq = working_seq

        return encoding_seq

    @property
    def one_hot_encode_sequence(self) -> Sequence[Sequence[int]]:
        """
        Encodes D x N where D is the size of the alphabet and N is the padding.

        Returns:
            A one hot encoded matrix representing the sequence.
        """

        encoded_sequence = self.alphabet.integer_encode(self._seq)
        encoded_len = len(encoded_sequence)
        letters_len = len(self.alphabet.letters_and_tokens)

        one_hot_seq = np.zeros((encoded_len, letters_len), dtype=np.int)
        one_hot_seq[np.arange(encoded_len), encoded_sequence] = 1

        return one_hot_seq.tolist()
