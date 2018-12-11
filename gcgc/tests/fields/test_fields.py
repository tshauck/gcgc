# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import unittest
from pathlib import Path

from gcgc.fields.categorical_field import FileMetaDataField


class TestFileMetaDataField(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_generate_encodings(self):

        paths = [Path("a"), Path("a.2"), Path("b")]
        expected_labels = {"a", "b"}

        def preprocess(p: Path) -> str:
            s = str(p)
            return s[0]

        ff = FileMetaDataField.from_paths("test", paths, preprocess)

        keys = set(ff.encoding_dict.keys())

        self.assertEqual(expected_labels, keys)
