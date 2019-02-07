# Rollouts

For some applications the sequence is not of direct interest, but rather some subsections of the
sequence. GCGC supports breaking up the sequence into some subsets then iterating over them, in what
is termed a rollout.

Right now, GCGC can rollout k-mers and sequence features that exist on the `Bio.SeqRecord.SeqRecord`.

## Rolling-out Kmers

Often sequences are broken up into "kmers", where "k" refers to a fixed length of the sequence. This
is somewhat analogous to ngrams in the NLP context.

For example, given an encoded sequence:

```python
from gcgc.rollouts import rollout_kmers

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

encoded_seq = SeqRecord(Seq("ATCGATCG"))
```

It is "rolled-out" with:

```python
rollout_iters = rollout_kmers(encoded_seq, kmer_length=4)
for rollout_kmer in rollout_iters:
    print(rollout_kmer)

# RolloutOutEncodedSeq(kmer=EncodedSeq('ATCG', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('A', IUPACUnambiguousDNAEncoding()))
# RolloutOutEncodedSeq(kmer=EncodedSeq('TCGA', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('T', IUPACUnambiguousDNAEncoding()))
# RolloutOutEncodedSeq(kmer=EncodedSeq('CGAT', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('C', IUPACUnambiguousDNAEncoding()))
# RolloutOutEncodedSeq(kmer=EncodedSeq('GATC', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('G', IUPACUnambiguousDNAEncoding()))
```

The sequence is of size 8 and the desired length is 4, so 4 sequences are produced. Also note that
optionally larger sizes of prior and next kmers can be included.

## SeqFeatures

SeqFeatures on SeqRecord are used to annotate things like genes or promoters. In addition to just
rolling out the features, an optional select function can be supplied to function that only will
yield the feature if it evaluate to True. The default is to select everything.

As an example, to create the SeqRecord with seq features.

```python
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


sr = SeqRecord(
    Seq("ATCG"),
    features=[
        SeqFeature(FeatureLocation(0, 4), type="gene"),
        SeqFeature(FeatureLocation(0, 4), type="promoter"),
    ],
)
```

Then roll it out only for genes.

```python
from gcgc.rollout import rollout_seq_features
rollout_features = rollout_seq_features(sr, lambda x: x.type == 'gene')

for rf in rollout_features:
    print(rf)

# Note the single output in light of two features.
# RolledOutEncodedSeqs(encoded_seq=EncodedSeq('ATCG', ExtendedIUPACProteinEncoding()), prior_encoded_seq=None, next_encoded_seq=None)

```
