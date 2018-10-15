# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from Bio import SeqIO


class FASTARecordGenerator(object):
    """
    """

    def __init__(self, directory_path, encoder):
        """
        """

        self.directory_path = directory_path
        self.encoder = encoder

    def __iter__(self):
        """
        """

        for file_path in self.directory_path.glob("*.fasta"):
            with open(file_path) as f:
                for seq_record in SeqIO.parse(f, "fasta"):
                    yield self.encoder.encode(seq_record)
