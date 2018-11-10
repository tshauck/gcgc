# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from typing import NamedTuple, Sequence

import tensorflow as tf
import numpy as np


def to_tensorflow_record(encoded_seq):
    """
    Convert the sequence to a tensorflow record.
    """

    integer_encoded = tf.train.Feature(
        int64_list=tf.train.Int64List(value=encoded_seq.integer_encoded)
    )

    alphabet_letters = tf.train.Feature(
        bytes_list=tf.train.BytesList(
            value=[m.encode("utf-8") for m in encoded_seq.alphabet.letters_and_tokens]
        )
    )

    feature_arg = {"integer_encoded": integer_encoded, "alphabet_letters": alphabet_letters}
    example = tf.train.Example(features=tf.train.Features(feature=feature_arg))

    return example


# TODO(trent): Fully functional EncodedSeq?
class ParsedEncodedSequence(NamedTuple):
    integer_encoded: np.ndarray
    alphabet_letters: np.ndarray


def from_tensorflow_example(example):
    features = {
        "integer_encoded": tf.FixedLenSequenceFeature((), tf.int64, allow_missing=True),
        "alphabet_letters": tf.FixedLenSequenceFeature((), tf.string, allow_missing=True),
    }
    parsed_features = tf.parse_single_example(example, features)
    return ParsedEncodedSequence(
        integer_encoded=parsed_features["integer_encoded"],
        alphabet_letters=parsed_features["alphabet_letters"],
    )
