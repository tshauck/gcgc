# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from pathlib import Path
from typing import Dict

import torch
from Bio.SeqRecord import SeqRecord

from gcgc.parser import SequenceParser


class TorchSequenceParser(SequenceParser):
    def parse_record(self, input_seq: SeqRecord, path: Path) -> Dict:
        """
        Convert the incoming SeqRecord to a dictionary of features.
        """

        parsed_features = super().parse_record(input_seq, path)
        parsed_features["seq_tensor"] = torch.LongTensor(parsed_features["seq_tensor"])

        parsed_features["id"] = input_seq.id

        if self.has_file_features:
            for ff in self.file_features:
                parsed_features[ff.name] = torch.tensor(parsed_features[ff.name])

        return parsed_features
