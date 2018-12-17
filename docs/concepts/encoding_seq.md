# Encoding Seq

The `EncodedSeq` object is similar to the BioPython `Seq` object, in that it primarily contains the
sequence letters and an associated Alphabet with the prime difference being the Alphabet is a GCGC
alphabet and not from BioPython.

## Creating an EncodedSeq

To create an `EncodedSeq` pass a sequence and an alphabet.

```python
from gcgc.encoded_seq import EncodedSeq
from gcgc.alphabet import IUPACUnambiguousDNAEncoding

es = EncodedSeq("ATCG", IUPACUnambiguousDNAEncoding())
# EncodedSeq('ATCG', IUPACUnambiguousDNAEncoding())
```

## Modifying an EncodedSeq

Once the object has been created, there are various ways to modify the underlying sequence.

```python
es.pad(pad_to=10)
# EncodedSeq('ATCG||||||', IUPACUnambiguousDNAEncoding())

es.encapsulate()
# EncodedSeq('>ATCG<', IUPACUnambiguousDNAEncoding())
```

The `EncodedSeq` object also supports chaining.

```python
es.encapsulate().conform(7)
# EncodedSeq('>ATCG<|', IUPACUnambiguousDNAEncoding())
```

## Integer Encodings

After the sequence has been modified, integer encodings are available as properties.

```python
es.one_hot_encoded
# [[0, 1, 0, 0, 0, 0, 0],
#  [0, 0, 1, 0, 0, 0, 0],
#  [0, 0, 0, 1, 0, 0, 0],
#  [1, 0, 0, 0, 0, 0, 0]]

es.integer_encoded
# [1, 2, 3, 0]
```

## Rolling-out Kmers

Often sequences are broken up into "kmers", where "k" refers to a fixed length of the sequence. This
is somewhat analogous to ngrams in the NLP context.

For example, given an encoded sequence:

```python
encoded_seq = EncodedSeq("ATCGATCG", IUPACUnambiguousDNAEncoding())
```

It is "rolled-out" with:

```python
rollout_iters = encoded_seq.rollout_kmers(kmer_length=4)
for rollout_kmer in rollout_iters:
    print(rollout_kmer)

# RolledOutKmer(kmer=EncodedSeq('ATCG', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('A', IUPACUnambiguousDNAEncoding()))
# RolledOutKmer(kmer=EncodedSeq('TCGA', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('T', IUPACUnambiguousDNAEncoding()))
# RolledOutKmer(kmer=EncodedSeq('CGAT', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('C', IUPACUnambiguousDNAEncoding()))
# RolledOutKmer(kmer=EncodedSeq('GATC', IUPACUnambiguousDNAEncoding()), prior_kmer=EncodedSeq('', IUPACUnambiguousDNAEncoding()), next_kmer=EncodedSeq('G', IUPACUnambiguousDNAEncoding()))
```

The sequence is of size 8 and the desired length is 4, so 4 sequences are produced. Also note that
optionally larger sizes of prior and next kmers can be included.
