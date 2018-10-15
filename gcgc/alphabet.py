# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from typing import Dict, Any, Sequence, NamedTuple
from abc import ABC

from Bio.Data import IUPACData
from Bio.SeqRecord import SeqRecord
import numpy as np


class EncodedSeqRecord(NamedTuple):
    """
    An inteface that holds information on encoded types.
    """

    one_hot_sequence: Sequence[Sequence[int]]
    sequence_id: str
    integer_encode: Sequence[int]


class EncodingAlphabet(ABC):
    """
    The Encoding Alphabet is meant to be a baseclass for other alphabets. This is similar to
    vocabularies often found in NLP packages, though the extent of the vocabulary is known upfront.
    """

    START = ">"
    END = "<"
    PADDING = "|"

    def __init__(self, letters, encapsulate: bool = True, padding_to: int = 500):
        """
        Create the EncodingAlphabet object.

        Args:
            letters: A list of tokens which constitue the alphabet.
            padding_to: When encoding a sequence, how long to pad a sequence to
                when encoding.
        """

        self.letters = letters
        self.letters_and_tokens = self.letters + self.START + self.END + self.PADDING
        self.encapsulate = encapsulate

        self.encoding_index = {
            letter: idx for idx, letter in enumerate(self.letters_and_tokens)
        }
        self.decoding_index = {
            idx: letter for letter, idx in self.encoding_index.items()
        }

        self.padding_to = padding_to

    def _seq(self, seq: str) -> str:
        """
        Return the sequence that is encapsulated in the tokens, if the associated flag is set.

        Args:
            seq: The alphabet of the sequence.

        """

        if self.encapsulate:
            return self.encapsulate_sequence(seq)
        else:
            return seq

    def encode_token(self, token: str) -> int:
        """
        Given a particular token, return the integer representation.

        Args:
            token: The token to convert.

        Returns:
            The integer representing the token.
        """

        return self.encoding_index[token]

    def decode_token(self, int_token: int) -> str:
        """
        Decode a token. This is the inverse of encode_token.

        Args:
            int_token: The integer representation of a token.

        Returns:
            The str which was encoded by the integer.
        """

        return self.decoding_index[int_token]

    def pad(self, seq) -> str:
        """
        Pad a sequence up to self.padding_to if it's shorter, otherwise return the sequence.

        Args:
            seq: The sequence to pad.

        Returns:
            The sequence padding up to padding_to characters.
        """

        working_seq = self._seq(seq)
        seq_len = len(working_seq)

        if seq_len < self.padding_to:
            extra_chars = self.padding_to - seq_len
            encoding_seq = working_seq + (self.PADDING * extra_chars)
        else:
            encoding_seq = working_seq

        return encoding_seq

    def integer_encode(self, seq: str) -> Sequence[int]:
        """
        Integer encode the sequence.

        Args:
            seq: The sequence to pad.

        Returns:
            The integer sequence representation of the sequence.
        """

        padded_seq = self.pad(seq)
        return np.array([self.encoding_index[s] for s in padded_seq])

    def one_hot_encode_sequence(self, seq: str) -> Sequence[Sequence[int]]:
        """
        Encodes D x N where D is the size of the alphabet and N is the padding.
        """

        # This is a bit jenky as even if we aren't encapselating the sequence
        # it still is used to determine size of the onehot encoding.

        encoded_sequence = self.integer_encode(seq)
        encoded_len = len(encoded_sequence)
        letters_len = len(self.letters_and_tokens)

        one_hot_seq = np.zeros((encoded_len, letters_len), dtype=np.int)
        one_hot_seq[np.arange(encoded_len), encoded_sequence] = 1

        return one_hot_seq.tolist()

    def encapsulate_sequence(self, seq: str) -> str:
        """
        Wrap a sequence with the start and end tokens.

        Args:
            seq: The sequence.

        Returns:
            The sequence wrapped in START and END.
        """
        return f"{self.START}{seq}{self.END}"

    def encode(self, seq: SeqRecord) -> EncodedSeqRecord:
        """
        Encode the sequence into an EncodedSeqRecord.

        Args:
            seq: The BioPython SeqRecord to encode.

        Returns:
            The SeqRecord encoded as an EncodedSeqRecord.
        """

        return EncodedSeqRecord(
            one_hot_sequence=self.one_hot_encode_sequence(seq.seq),
            sequence_id=seq.id,
            integer_encoded=self.integer_encode(seq.seq),
        )


class IUPACProteinExtendedAlphabet(EncodingAlphabet):
    def __init__(self, encapsulate_sequence: bool = True, padding_to: int = 50):
        super(IUPACProteinExtendedAlphabet, self).__init__(
            IUPACData.extended_protein_letters, encapsulate_sequence, padding_to
        )


class IUPACDnaExtendedAlphabet(EncodingAlphabet):
    def __init__(self, encapsulate_sequence: bool = True, padding_to: int = 50):
        super(IUPACDnaExtendedAlphabet, self).__init__(
            IUPACData.extended_dna_letters, encapsulate_sequence, padding_to
        )


class UnambiguousDnaExtendedAlphabet(EncodingAlphabet):
    def __init__(self, encapsulate_sequence: bool = True, padding_to: int = 50):
        super(UnambiguousDnaExtendedAlphabet, self).__init__(
            IUPACData.unambiguous_dna_letters, padding_to=padding_to
        )
