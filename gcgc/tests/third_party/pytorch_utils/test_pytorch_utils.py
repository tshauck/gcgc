# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import tempfile
import pathlib

import pytest
from torch.utils.data import DataLoader

from gcgc.alphabet.iupac import IUPACProteinEncoding
from gcgc.ml.pytorch_utils.data import GenomicDataset
from gcgc.parser import SequenceParser
from gcgc.tests.fixtures import ECOLI_PATH

SP = SequenceParser()
POST_PAD_LENGTH = 1000


@pytest.yield_fixture
def genomic_dataset():
    glob = sorted(ECOLI_PATH.glob("*.fasta"))

    pe = IUPACProteinEncoding()

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = pathlib.Path(tmpdir) / "test.db"
        yield GenomicDataset.from_paths(glob, SP, "fasta", pe, str(db_path), POST_PAD_LENGTH)


def test_index_multiple_files(genomic_dataset):

    assert len(genomic_dataset) == 25

    test_sequences = {
        0: "sp|C5A0C3|ARGE_ECOBW",
        1: "sp|P23908|ARGE_ECOLI",
        24: "sp|B1X825|LFTR_ECODH",
    }

    for idx, expected_id in test_sequences.items():
        actual_record = genomic_dataset[idx]
        assert actual_record["id"] == expected_id
        assert len(actual_record["seq_tensor"]) == POST_PAD_LENGTH


@pytest.mark.skip("Skipping this until it can be fixed, issue with random access to files.")
def test_index_num_workers(genomic_dataset):
    dataloader = DataLoader(genomic_dataset, num_workers=4)
    for i in range(10):
        for _ in dataloader:
            pass
