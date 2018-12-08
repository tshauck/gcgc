# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import torch
from Bio.SeqRecord import SeqRecord

from gcgc.encoded_seq import EncodedSeq
from gcgc.exceptions import EncodedSeqLengthParserException
from gcgc.fields import FileMetaDataField


@dataclass
class EncodedSeqLengthParser:
    pad_to: Optional[int] = None
    conform_to: Optional[int] = None

    def parse_encoded_seq_record(self, encoded_seq: EncodedSeq) -> EncodedSeq:
        num_non_none = sum(1 for x in [self.pad_to, self.conform_to] if x is not None)
        if num_non_none > 1:
            raise EncodedSeqLengthParserException("One one of pad_to or conform_to an be set.")

        if self.pad_to:
            return encoded_seq.pad(pad_to=self.pad_to)

        if self.conform_to:
            return encoded_seq.conform(conform_to=self.conform_to)

        raise RuntimeError(
            "Could not parse seq length for some reason, check parser configuration."
            f"conform_to={self.conform_to} pad_to={self.pad_to}"
        )


@dataclass
class SequenceParser:
    """
    A configurable sequence parser.
    """

    encapsulate: bool = True
    seq_length_parser: Optional[EncodedSeqLengthParser] = None
    file_features: Optional[List[FileMetaDataField]] = None

    def parse_record(self, input_seq: SeqRecord, path: Path) -> Dict:
        """
        Convert the incoming SeqRecord to a dictionary of features.
        """

        parsed_features: Dict = {}

        if self.encapsulate:
            es = EncodedSeq.from_seq(input_seq.seq).encapsulate()

        if self.seq_length_parser:
            es = self.seq_length_parser.parse_encoded_seq_record(es)

        parsed_features["seq_tensor"] = torch.LongTensor(es.integer_encoded)
        parsed_features["id"] = input_seq.id

        parsed_features.update(self._generate_file_features(path))

        return parsed_features

    @property
    def has_file_features(self) -> bool:
        return bool(self.file_features)

    def _generate_file_features(self, path):
        file_features = {}

        if self.has_file_features:
            for ff in self.file_features:
                file_features[ff.name] = ff.torch_tensor(path)

        return file_features
