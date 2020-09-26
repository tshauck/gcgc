# GCGC

> GCGC is a tool for feature processing on Biological Sequences.

[![](https://github.com/tshauck/gcgc/workflows/Run%20Tests%20and%20Lint/badge.svg)](https://github.com/tshauck/gcgc/actions?query=workflow%3A%22Run+Tests+and+Lint%22)
[![](https://img.shields.io/pypi/v/gcgc.svg)](https://pypi.python.org/pypi/gcgc)
[![code style black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Installation

GCGC is primarily intended to be used as part of a larger workflow inside
Python.

To install via pip:

```sh
$ pip install gcgc
```

If you'd like to use code that helps gcgc's tokenizers integrate with common
third party libraries, either install those packages separately, or use gcgc's
extras.

```sh
$ pip install 'gcgc[pytorch,hf]'
```

## Documentation

The GCGC documentation is at [gcgc.trenthauck.com](http://gcgc.trenthauck.com),
please see it for examples.

### Quick Start

The easiest way to get started is to import the kmer tokenizer, configure it,
then start tokenizing.

```python
from gcgc import KmerTokenizer

kmer_tokenizer = KmerTokenizer(alphabet="unambiguous_dna")
encoded = kmer_tokenizer.encode("ATCG")
print(encoded)
```

sample output:

```
[1, 6, 7, 8, 5, 2]
```

This output includes the "bos" token, the "eos" token, and the four nucleotide
tokens in between.

You can go the other way and convert the integers to strings.

```python
from gcgc import KmerTokenizer

kmer_tokenizer = KmerTokenizer(alphabet="unambiguous_dna")
decoded = kmer_tokenizer.decode(kmer_tokenizer.encode("ATCG"))
print(decoded)
```

sample output:

```
['>', 'A', 'T', 'C', 'G', '<']
```

There's also the vocab for the kmer tokenizer.

```python
from gcgc import KmerTokenizer

kmer_tokenizer = KmerTokenizer(alphabet="unambiguous_dna")
print(kmer_tokenizer.vocab.stoi)
```

sample output:

```
{'|': 0, '>': 1, '<': 2, '#': 3, '?': 4, 'G': 5, 'A': 6, 'T': 7, 'C': 8}
```

### Transformers

The [Transformers](https://huggingface.co/transformers/) library has an
idea of a tokenizer that is used for various modeling tasks.

To make it easier to use the Transformers library on biological sequences, gcgc
has a Transformers compatible tokenizer that can be created from the
KmerTokenizer.

```python
from gcgc import KmerTokenizer
from gcgc.third_party.hf import GCGCTransformersTokenizer

kmer_tokenizer = KmerTokenizer(
  kmer_length=2, kmer_stride=2, alphabet="unambiguous_dna"
)

tt = GCGCTransformersTokenizer.from_kmer_tokenizer(
    kmer_tokenizer
)

batch = tt.batch_encode_plus(["ATGC", "GCGC"])
print(batch["input_ids"])
```

sample output:

```
[[1, 11, 8, 2], [1, 8, 8, 2]]
```
