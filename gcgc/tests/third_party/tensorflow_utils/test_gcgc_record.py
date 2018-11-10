# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest

from Bio import SeqIO

import tensorflow as tf

from gcgc.encoded_seq import EncodedSeq
from gcgc.alphabet.iupac import ExtendedIUPACDNAEncoding
from gcgc.third_party.tensorflow_utils import record as gcgc_record
from gcgc.tests.fixtures import P53_HUMAN


class TestTFRecordSerialization(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_tf_records(self):

        with open(P53_HUMAN) as f:
            sr = SeqIO.read(f, "fasta")

        encoded_seq = EncodedSeq.from_seq(sr.seq)
        filenames = [str(P53_HUMAN.with_suffix(".tf_records"))]

        dataset = tf.data.TFRecordDataset(filenames)
        dataset = dataset.map(gcgc_record.from_tensorflow_example)
        dataset_iterator = dataset.make_one_shot_iterator()
        next_item = dataset_iterator.get_next()

        session = tf.Session()
        seq = session.run(next_item)

        self.assertEqual(encoded_seq.integer_encoded, seq.integer_encoded.tolist())
