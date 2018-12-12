# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from typing import Dict

import torch

from gcgc.parser import SequenceParser
from gcgc.parser.gcgc_record import GCGCRecord


class TorchSequenceParser(SequenceParser):
    def parse_record(self, gcgc_record: GCGCRecord) -> Dict:
        """
        Convert the incoming SeqRecord to a dictionary of features.
        """

        parsed_features = super().parse_record(gcgc_record)
        parsed_features["seq_tensor"] = torch.LongTensor(parsed_features["seq_tensor"])

        if self.has_file_features:
            for ff in self.file_features:
                parsed_features[ff.name] = torch.tensor(parsed_features[ff.name])

        return parsed_features
