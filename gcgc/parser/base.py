# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A Parser that converts GCGCRecords into data suitible for ML training."""

from typing import Any, Dict, List, Optional

from gcgc.encoded_seq import EncodedSeq
from gcgc.exceptions import EncodedSeqLengthParserException
from gcgc.fields import AnnotationField, DescriptionField, FileMetaDataField
from gcgc.parser.gcgc_record import GCGCRecord


class EncodedSeqLengthParser:
    """A parser to turn an input EncodedSeq into one of the specified length."""

    def __init__(self, pad_to: Optional[int] = None, conform_to: Optional[int] = None) -> None:
        """Create a parser to work with sequence length."""

        self.pad_to = pad_to
        self.conform_to = conform_to

    def parse_encoded_seq_record(self, encoded_seq: EncodedSeq) -> EncodedSeq:
        """Return a new EncodedSeq object at the proper length given the specification."""

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


class SequenceParser:
    """A configurable sequence parser."""

    def __init__(
        self,
        encapsulate: bool = True,
        seq_length_parser: Optional[EncodedSeqLengthParser] = None,
        file_features: Optional[List[FileMetaDataField]] = None,
        annotation_features: Optional[List[AnnotationField]] = None,
        description_features: Optional[List[DescriptionField]] = None,
        sequence_offset: Optional[int] = None,
    ) -> None:
        """Create the SequenceParser object."""

        self.encapsulate = encapsulate
        self.seq_length_parser = seq_length_parser

        self.file_features = file_features if file_features is not None else []
        self.annotation_features = annotation_features if annotation_features is not None else []
        self.description_features = description_features if description_features is not None else []

        self.sequence_offset = sequence_offset

    def _preprocess_record(self, es: EncodedSeq):
        if self.encapsulate:
            es = es.encapsulate()

        if self.seq_length_parser:
            es = self.seq_length_parser.parse_encoded_seq_record(es)

        return es

    def parse_record(self, gcgc_record: GCGCRecord) -> Dict:
        """Convert the incoming GCGCRecord to a dictionary of features."""

        es = gcgc_record.encoded_seq
        processed_seq = self._preprocess_record(es)

        parsed_features: Dict[str, Any] = {}

        parsed_features["seq_tensor"] = processed_seq.integer_encoded
        parsed_features["seq_tensor_one_hot"] = processed_seq.one_hot_encoded

        if self.has_offset:
            offset_seq = processed_seq.shift(self.sequence_offset)

            parsed_features["offset_seq_tensor"] = offset_seq.integer_encoded
            parsed_features["offset_seq_tensor_one_hot"] = offset_seq.one_hot_encoded

        parsed_features["id"] = gcgc_record.seq_record.id
        parsed_features.update(self._generate_file_features(gcgc_record.path))
        parsed_features.update(self._generate_annotation_features(gcgc_record.seq_record))
        parsed_features.update(self._generate_description_features(gcgc_record.seq_record))

        return parsed_features

    @property
    def has_description_features(self) -> bool:
        """Return True if this parser has description features."""
        return bool(self.description_features)

    @property
    def has_offset(self) -> bool:
        """Return True if this object has an offset."""
        return self.sequence_offset is not None

    @property
    def has_file_features(self) -> bool:
        """Return True if this parser has file features."""
        return bool(self.file_features)

    def _generate_file_features(self, path):
        file_features = {}

        if self.has_file_features:
            for file_feature in self.file_features:
                file_features[file_feature.name] = file_feature.encode(path)

        return file_features

    @property
    def has_annotation_features(self) -> bool:
        """Return True if this parser has annotation features."""
        return bool(self.annotation_features)

    def _generate_annotation_features(self, seq_record):
        seq_records_features = {}

        if self.has_annotation_features:
            for ff in self.annotation_features:
                seq_records_features[ff.name] = ff.encode(seq_record.annotations)

        return seq_records_features

    def _generate_description_features(self, seq_record):
        description_features = {}

        if self.has_description_features:
            for df in self.description_features:
                description_features[df.name] = df.encode(seq_record.description)

        return description_features
