# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Base rollout module."""

from typing import Callable, Generator, Optional

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from gcgc.encoded_seq import EncodedSeq


class RolledOutEncodedSeqs:
    """Holds information on the rollout context."""

    def __init__(
        self,
        encoded_seq: EncodedSeq,
        prior_encoded_seq: Optional[EncodedSeq] = None,
        next_encoded_seq: Optional[EncodedSeq] = None,
    ):
        """Create the rolled out encoded sequence object."""
        self.encoded_seq = encoded_seq
        self.prior_encoded_seq = prior_encoded_seq
        self.next_encoded_seq = next_encoded_seq

    def apply(self, func: Callable[[EncodedSeq], EncodedSeq]) -> "RolledOutEncodedSeqs":
        """Apply the callable func to each sequence and return a copy of the object."""

        def if_null(es: Optional[EncodedSeq] = None):
            """Handle the None case of applied to func."""
            return func(es) if es is not None else None

        return RolledOutEncodedSeqs(
            if_null(self.encoded_seq),
            if_null(self.prior_encoded_seq),
            if_null(self.next_encoded_seq),
        )


def rollout_kmers(
    sr: SeqRecord,
    kmer_length: int = 5,
    prior_length: int = 0,
    next_kmer_length: int = 1,
    window: int = 1,
    func: Optional[Callable[[EncodedSeq], EncodedSeq]] = None,
) -> Generator[RolledOutEncodedSeqs, None, None]:
    """Rollout kmers of length k from the sequence as RolledOutEncodedSeqs."""

    es = EncodedSeq.from_seq(sr.seq)

    seq_length = len(es)

    rollout_start = prior_length
    rollout_to = seq_length - kmer_length

    for i in range(rollout_start, rollout_to, window):
        prior_kmer = es[(i - prior_length) : i]
        kmer = es[i : i + kmer_length]
        next_kmer = es[i + kmer_length : i + kmer_length + next_kmer_length].conform(
            next_kmer_length
        )

        rollout_kmer = RolledOutEncodedSeqs(kmer, prior_kmer, next_kmer)
        if func is not None:
            rollout_kmer = rollout_kmer.apply(func)

        yield rollout_kmer


def default_select_func(sf: SeqFeature) -> bool:
    """Return True by default."""
    return True


def rollout_seq_features(
    sr: SeqRecord, select_func: Callable[[SeqFeature], bool] = default_select_func
) -> Generator[RolledOutEncodedSeqs, None, None]:
    """Rollout sequence features as RolledOutEncodedSeqs."""

    for sf in sr.features:
        if select_func(sf):
            sf_sr = sf.extract(sr)
            sf_es = EncodedSeq.from_seq(sf_sr.seq)

            yield RolledOutEncodedSeqs(sf_es)
