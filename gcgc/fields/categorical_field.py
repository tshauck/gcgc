# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List

from gcgc.fields.field import Field


@dataclass
class LabelField(Field):
    """
    A type of Field that can work with label values -- roundtripping them to an integer encoding.
    """

    encoding_dict: Dict[str, int]
    decoding_dict: Dict[int, str]

    def encode(self, label):
        """
        Look up the label in the encoding dict and return it.
        """

        return self.encoding_dict[label]

    def decode(self, label_int):
        """
        Look up the label's integer from the decoding dict and return it.
        """

        return self.decoding_dict[label_int]

    @classmethod
    def from_vocabulary(cls, name: str, vocab: List[str]):
        encoding_dict = {}
        decoding_dict = {}

        for i, s in enumerate(sorted(set(vocab))):
            encoding_dict[s] = i
            decoding_dict[i] = s

        return LabelField(name, encoding_dict, decoding_dict)


def default_preprocess(p: Path) -> str:
    return str(p)


@dataclass
class FileMetaDataField(LabelField):
    """
    """

    preprocess: Callable[[Path], str] = default_preprocess

    def encode(self, file: Path) -> int:
        """
        Preprocess the file path, then encode the resultant label.
        """
        label = self.preprocess(file)
        return super().encode(label)

    @classmethod
    def from_paths(
        cls, name: str, paths: List[Path], preprocess: Callable[[Path], str] = default_preprocess
    ):
        """
        Given a set of exemplar paths, create the (d-)encoding dict and return the field.
        """

        str_vocab = [preprocess(p) for p in paths]
        label_field = super().from_vocabulary(name, str_vocab)

        return cls(name, label_field.encoding_dict, label_field.decoding_dict, preprocess)


@dataclass
class AnnotationField(LabelField):
    """
    """

    preprocess: Callable[[Dict], str]

    def encode(self, file: Path) -> int:
        """
        Preprocess the file path, then encode the resultant label.
        """
        label = self.preprocess(file)
        return super().encode(label)

    @classmethod
    def from_annotations(
        cls, name: str, annotations: List[Dict], preprocess: Callable[[Dict], str]
    ):
        """
        Given a set of exemplar annotations, create the (d-)encoding dict and return the field.
        """

        str_vocab = [preprocess(a) for a in annotations]
        label_field = super().from_vocabulary(name, str_vocab)

        return cls(name, label_field.encoding_dict, label_field.decoding_dict, preprocess)
