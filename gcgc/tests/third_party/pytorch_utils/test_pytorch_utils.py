# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest

from gcgc.alphabet.iupac import IUPACProteinEncoding
from gcgc.parser import SequenceParser
from gcgc.tests.fixtures import ECOLI_PATH, P53_HUMAN
from gcgc.third_party.pytorch_utils.data import GenomicDataset


class TestPyTorchUtils(unittest.TestCase):
    def setUp(self) -> None:
        self.sp = SequenceParser()

    def test_load_dataset(self):
        def yielder():
            yield P53_HUMAN

        test_dataset = GenomicDataset.from_paths(yielder(), self.sp, "fasta")
        self.assertEqual(len(test_dataset), 1)

    def test_index_multiple_files(self):

        glob = ECOLI_PATH.glob("*.fasta")

        pe = IUPACProteinEncoding()
        test_dataset = GenomicDataset.from_paths(glob, self.sp, "fasta", pe)
        self.assertEqual(len(test_dataset), 25)

        test_sequences = {
            0: "sp|C5A0C3|ARGE_ECOBW",
            1: "sp|P23908|ARGE_ECOLI",
            24: "sp|B1X825|LFTR_ECODH",
        }

        for idx, expected_id in test_sequences.items():
            actual_record = test_dataset[idx]
            self.assertEqual(actual_record["id"], expected_id)
