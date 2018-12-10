# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import tempfile
import unittest
from pathlib import Path

from gcgc.datasets import dataset


class TestDataset(unittest.TestCase):
    def test_uniprot_dataset(self):

        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdir_path = Path(tmpdirname)
            taxon_dataset = dataset.TaxonDataset([83333], tmpdir_path)
            taxon_dataset.download_files()

            files = list(tmpdir_path.glob("*.fasta"))

            self.assertEqual(len(files), 1)
            self.assertTrue(all(f.stat().st_size > 0 for f in files))
