# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Test a keras model works with the tokenization."""

import tensorflow as tf

import pathlib
import csv

from gcgc import SequenceTokenizer, SequenceTokenizerSpec

HERE = pathlib.Path(__file__).parent


CLASSES = {"N": 2, "EI": 0, "IE": 1}


def get_dataset(tokenizer):
    # https://archive.ics.uci.edu/ml/machine-learning-databases/molecular-biology/splice-junction-gene-sequences/
    split_tsv = HERE / "./split.tsv"

    def yield_records():
        with split_tsv.open("r") as tsv_handler:
            rows = csv.DictReader(tsv_handler, fieldnames=["label", "sequence"], delimiter="\t")
            for row in rows:
                yield {"encoded_seq": tokenizer.encode(row["sequence"])}, CLASSES[row["label"]]

    return tf.data.Dataset.from_generator(
        yield_records,
        ({"encoded_seq": tf.int32}, tf.int32),
        (
            {"encoded_seq": tf.TensorShape([tokenizer.tokenizer_spec.max_length])},
            tf.TensorShape([]),
        ),
    ).batch(64)


class CNNClassifier(tf.keras.Model):
    def __init__(self, tokenizer, *args, **kwargs):
        super(CNNClassifier, self).__init__(*args, **kwargs)
        self.tokenizer = tokenizer
        self.embedding = tf.keras.layers.Embedding(
            len(self.tokenizer.tokenizer_spec.vocabulary), 50
        )
        self.cnn = tf.keras.layers.Conv1D(50, 10)
        self.cnn2 = tf.keras.layers.Conv1D(50, 10)

        self.flatten = tf.keras.layers.Flatten()
        self.out = tf.keras.layers.Dense(500)
        self.out2 = tf.keras.layers.Dense(3)

    def call(self, inputs, **kwargs):
        seq = inputs["encoded_seq"]

        embedded = self.embedding(seq)
        cnn = self.cnn2(self.cnn(embedded))
        flattened = self.flatten(cnn)
        return self.out2(tf.nn.relu(self.out(flattened)))


def test_fit_model():
    """Test we can fit a keras model."""
    max_length = 60

    spec = SequenceTokenizerSpec(max_length, "ATCGNDRS")
    tokenizer = SequenceTokenizer(spec)

    optimizer = tf.keras.optimizers.Adam(learning_rate=3e-6, epsilon=1e-08, clipnorm=1.0)
    loss = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
    metric = tf.keras.metrics.SparseCategoricalAccuracy("accuracy")

    dataset = get_dataset(tokenizer)
    model = CNNClassifier(tokenizer)
    model.compile(loss=loss, optimizer=optimizer, metrics=[metric])
    model.fit(dataset, epochs=1)
