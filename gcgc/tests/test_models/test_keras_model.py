# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Test a keras model works with the tokenization."""

import pathlib
import csv

import tensorflow as tf

from gcgc import SequenceTokenizer, SequenceTokenizerSpec

HERE = pathlib.Path(__file__).parent


CLASSES = {"N": 2, "EI": 0, "IE": 1}


def load_and_tokenize_dataset(tokenizer: SequenceTokenizer):
    """Load the dataset and create a tf.data.Dataset from the underlying records.

    Data from:
        https://archive.ics.uci.edu/ml/machine-learning-databases/molecular-biology/splice-junction-gene-sequences/

    Args:
        tokenizer: A tokenizer that converts the input inputs into a format suitable for keras.

    Returns:
        A tensorflow dataset suitable for model.fit.

    """
    split_tsv = HERE / "./split.tsv"

    def yield_records():
        """Yield encoded records."""
        with split_tsv.open("r") as tsv_handler:
            rows = csv.DictReader(tsv_handler, fieldnames=["label", "sequence"], delimiter="\t")
            for row in rows:
                yield {"encoded_seq": tokenizer.encode(row["sequence"])}, CLASSES[row["label"]]

    return (
        tf.data.Dataset.from_generator(
            yield_records,
            ({"encoded_seq": tf.int32}, tf.int32),
            (
                {"encoded_seq": tf.TensorShape([tokenizer.tokenizer_spec.max_length])},
                tf.TensorShape([]),
            ),
        )
        .shuffle(buffer_size=1000)
        .batch(64)
    )


class CNNClassifier(tf.keras.Model):
    """An object that runs a CNN over the sequences."""

    def __init__(self, tokenizer, *args, **kwargs):
        """Init the class."""
        super(CNNClassifier, self).__init__(*args, **kwargs)

        self.tokenizer = tokenizer

        self.submodel = tf.keras.models.Sequential(
            [
                tf.keras.layers.InputLayer(input_shape=(60,)),
                tf.keras.layers.Embedding(len(self.tokenizer.tokenizer_spec.vocabulary), 50),
                tf.keras.layers.Conv1D(12, 10),
                tf.keras.layers.Flatten(),
                tf.keras.layers.ReLU(),
                tf.keras.layers.Dense(500),
                tf.keras.layers.ReLU(),
                tf.keras.layers.Dense(3),
            ]
        )

    def call(self, inputs):
        """Implment the call method of the model that runs the forward pass."""
        seq = inputs["encoded_seq"]
        return self.submodel(seq)


def test_fit_model():
    """Test we can fit a keras model."""
    max_length = 60

    spec = SequenceTokenizerSpec("ATCGNDRS", max_length=max_length)
    tokenizer = SequenceTokenizer(spec)

    optimizer = tf.keras.optimizers.Adam()
    loss = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
    metric = tf.keras.metrics.SparseCategoricalAccuracy("accuracy")

    dataset = load_and_tokenize_dataset(tokenizer)
    model = CNNClassifier(tokenizer)
    model.compile(loss=loss, optimizer=optimizer, metrics=[metric])
    model.fit(dataset, epochs=25, verbose=2)
