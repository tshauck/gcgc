# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""PyTorch specific parser."""

from typing import Dict
from typing import Optional

import torch

from gcgc.parser import SequenceParser
from gcgc.parser.gcgc_record import GCGCRecord


class TorchSequenceParser(SequenceParser):
    """A PyTorch Sequence Parser."""

    def parse_record(self, gcgc_record: GCGCRecord, parsed_seq_len: Optional[int] = None) -> Dict:
        """Convert the incoming SeqRecord to a dictionary of features."""

        parsed_features = super().parse_record(gcgc_record, parsed_seq_len)
        parsed_features["seq_tensor"] = torch.LongTensor(parsed_features["seq_tensor"])
        parsed_features["seq_tensor_masked"] = torch.LongTensor(
            parsed_features["seq_tensor_masked"]
        )
        parsed_features["seq_len"] = torch.tensor(parsed_features["seq_len"])

        if self.has_offset:
            parsed_features["offset_seq_tensor"] = torch.LongTensor(
                parsed_features["offset_seq_tensor"]
            )

        if self.has_file_features:
            for file_feature in self.file_features:
                parsed_features[file_feature.name] = torch.tensor(
                    parsed_features[file_feature.name]
                )

        return parsed_features
