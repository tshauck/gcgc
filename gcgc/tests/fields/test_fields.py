# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from pathlib import Path
from typing import Dict
import unittest

from gcgc.fields.categorical_field import AnnotationField, FileMetaDataField


class TestFields(unittest.TestCase):
    def test_file_metadata_generate_encodings(self):

        paths = [Path("a"), Path("a.2"), Path("b")]
        expected_labels = {"a", "b"}

        def preprocess(p: Path) -> str:
            s = str(p)
            return s[0]

        ff = FileMetaDataField.from_paths("test", paths, preprocess)

        keys = set(ff.encoding_dict.keys())

        self.assertEqual(expected_labels, keys)

    def test_annotation_fields(self):

        annotations = [{"a": "a"}, {"a": "b"}]
        expected_labels = {"a", "b"}

        def preprocess(a: Dict) -> str:
            return a["a"]

        ff = AnnotationField.from_annotations("a", annotations, preprocess)

        keys = set(ff.encoding_dict.keys())

        self.assertEqual(expected_labels, keys)
