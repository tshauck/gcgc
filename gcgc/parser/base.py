# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from dataclasses import dataclass
from typing import Dict, List, Optional

from gcgc.encoded_seq import EncodedSeq
from gcgc.exceptions import EncodedSeqLengthParserException
from gcgc.fields import AnnotationField, FileMetaDataField
from gcgc.parser.gcgc_record import GCGCRecord


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
    annotation_features: Optional[List[AnnotationField]] = None

    def parse_record(self, gcgc_record: GCGCRecord) -> Dict:
        """
        Convert the incoming SeqRecord to a dictionary of features.
        """

        parsed_features: Dict = {}

        if self.encapsulate:
            es = EncodedSeq.from_seq(gcgc_record.seq_record.seq).encapsulate()

        if self.seq_length_parser:
            es = self.seq_length_parser.parse_encoded_seq_record(es)

        parsed_features["id"] = gcgc_record.seq_record.id
        parsed_features["seq_tensor"] = es.integer_encoded
        parsed_features.update(self._generate_file_features(gcgc_record.path))
        parsed_features.update(self._generate_annotation_features(gcgc_record.seq_record))

        return parsed_features

    @property
    def has_file_features(self) -> bool:
        return bool(self.file_features)

    def _generate_file_features(self, path):
        file_features = {}

        if self.has_file_features:
            for ff in self.file_features:
                file_features[ff.name] = ff.encode(path)

        return file_features

    @property
    def has_annotation_features(self) -> bool:
        return bool(self.annotation_features)

    def _generate_annotation_features(self, seq_record):
        seq_records_features = {}

        if self.has_annotation_features:
            for ff in self.annotation_features:
                seq_records_features[ff.name] = ff.encode(seq_record.annotations)

        return seq_records_features
