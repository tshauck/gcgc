# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

from typing import Callable, Dict, List, NamedTuple, Tuple
import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gcgc.alphabet.iupac import ExtendedIUPACDNAEncoding
from gcgc.rollout import rollout_kmers, rollout_seq_features


class TestRollOuts(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_seq_feature_rollout(self):
        class SeqFeatureTestSet(NamedTuple):
            name: str
            sr: SeqRecord
            expected_seqs: List[str]
            select_func: Callable[[SeqRecord], bool]

        test_table = [
            SeqFeatureTestSet(
                name="Test Full Seq",
                sr=SeqRecord(
                    Seq("ATCG", ExtendedIUPACDNAEncoding()),
                    features=[SeqFeature(FeatureLocation(0, 4))],
                ),
                expected_seqs=["ATCG"],
                select_func=lambda sr: True,
            ),
            SeqFeatureTestSet(
                name="Test Full Seq",
                sr=SeqRecord(
                    Seq("ATCG", ExtendedIUPACDNAEncoding()),
                    features=[
                        SeqFeature(FeatureLocation(0, 4), type="gene"),
                        SeqFeature(FeatureLocation(0, 4), type="promoter"),
                    ],
                ),
                expected_seqs=["ATCG"],
                select_func=lambda sr: sr.type == "gene",
            ),
        ]

        for test_set in test_table:
            rollout_features = rollout_seq_features(test_set.sr, test_set.select_func)

            for i, rf in enumerate(rollout_features):
                self.assertEqual(str(rf.encoded_seq), test_set.expected_seqs[i])

    def test_kmer_rollout(self):
        class KMerTestSet(NamedTuple):
            name: str
            sr: SeqRecord
            rollout_options: Dict[str, int]
            expected_seqs: List[Tuple[str, str, str]]

        test_table = [
            KMerTestSet(
                name="Test length of kmer",
                sr=SeqRecord(Seq("ATCGATCG", ExtendedIUPACDNAEncoding())),
                rollout_options={"kmer_length": 4},
                expected_seqs=[
                    ("ATCG", "", "A"),
                    ("TCGA", "", "T"),
                    ("CGAT", "", "C"),
                    ("GATC", "", "G"),
                ],
            ),
            KMerTestSet(
                name="Test length and prior length.",
                sr=SeqRecord(Seq("ATCGATCG", ExtendedIUPACDNAEncoding())),
                rollout_options={"kmer_length": 4, "prior_length": 2},
                expected_seqs=[("CGAT", "AT", "C"), ("GATC", "TC", "G")],
            ),
            KMerTestSet(
                name="Test Window",
                sr=SeqRecord(Seq("ATCGATCGATCG", ExtendedIUPACDNAEncoding())),
                rollout_options={"kmer_length": 4, "prior_length": 2, "window": 2},
                expected_seqs=[("CGAT", "AT", "C"), ("ATCG", "CG", "A"), ("CGAT", "AT", "C")],
            ),
            KMerTestSet(
                name="Test Window",
                sr=SeqRecord(Seq("ATCGATC", ExtendedIUPACDNAEncoding())),
                rollout_options={"kmer_length": 4, "next_kmer_length": 2},
                expected_seqs=[("ATCG", "", "AT"), ("TCGA", "", "TC"), ("CGAT", "", "C|")],
            ),
        ]

        for test_set in test_table:
            rollout_iters = rollout_kmers(test_set.sr, **test_set.rollout_options)
            for i, rk in enumerate(rollout_iters):

                expected_encoded_seq, expected_prior, expected_next = test_set.expected_seqs[i]

                self.assertEqual(str(expected_encoded_seq), rk.encoded_seq)
                self.assertEqual(str(expected_prior), rk.prior_encoded_seq)
                self.assertEqual(str(expected_next), rk.next_encoded_seq)
