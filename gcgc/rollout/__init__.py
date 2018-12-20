# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from typing import NamedTuple, Optional, Generator, Callable

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from gcgc.encoded_seq import EncodedSeq


class RolledOutEncodedSeqs(NamedTuple):
    encoded_seq: EncodedSeq
    prior_encoded_seq: Optional[EncodedSeq] = None
    next_encoded_seq: Optional[EncodedSeq] = None


def rollout_kmers(
    sr: SeqRecord,
    kmer_length: int = 5,
    prior_length: int = 0,
    next_kmer_length: int = 1,
    window: int = 1,
) -> Generator[RolledOutEncodedSeqs, None, None]:
    """Rollout kmers of length k from the sequence as RolledOutEncodedSeqs."""

    es = EncodedSeq.from_seq(sr.seq)

    seq_length = len(es)

    rollout_start = prior_length
    rollout_to = seq_length - kmer_length

    for i in range(rollout_start, rollout_to, window):
        prior_kmer = es[(i - prior_length) : i]
        kmer = es[i : i + kmer_length]
        next_kmer = es[i + kmer_length : i + kmer_length + next_kmer_length]

        rollout_kmer = RolledOutEncodedSeqs(kmer, prior_kmer, next_kmer)

        yield rollout_kmer


def default_select_func(sf: SeqFeature) -> bool:
    """Default func is to pass through."""
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
