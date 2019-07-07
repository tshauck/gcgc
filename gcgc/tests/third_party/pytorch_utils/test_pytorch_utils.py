# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import tempfile
import pathlib

from gcgc.alphabet.iupac import IUPACProteinEncoding
from gcgc.ml.pytorch_utils.data import GenomicDataset
from gcgc.parser import SequenceParser
from gcgc.tests.fixtures import ECOLI_PATH
from gcgc.tests.fixtures import P53_HUMAN

SP = SequenceParser()


def test_load_dataset():
    def yielder():
        yield P53_HUMAN

    test_dataset = GenomicDataset.from_paths(yielder(), SP, "fasta")
    assert len(test_dataset) == 1


def test_index_multiple_files():

    glob = sorted(ECOLI_PATH.glob("*.fasta"))

    pe = IUPACProteinEncoding()

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = pathlib.Path(tmpdir) / "test.db"
        test_dataset = GenomicDataset.from_paths(glob, SP, "fasta", pe, str(db_path))

        assert len(test_dataset) == 25

        test_sequences = {
            0: "sp|C5A0C3|ARGE_ECOBW",
            1: "sp|P23908|ARGE_ECOLI",
            24: "sp|B1X825|LFTR_ECODH",
        }

        for idx, expected_id in test_sequences.items():
            actual_record = test_dataset[idx]
            assert actual_record["id"] == expected_id
