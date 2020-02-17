# GCGC

GCGC is a package written in Python for pre-processing for biological sequences. Think of it like a
Natural Language Processing pre-processing toolkit with design choices oriented towards the
differences the sequences found in natural language vs biology.

At its core is a tokenizer specifically for biological sequences that can be
used as part of PyTorch or TensorFlow workflows.

## Installation

Install GCGC via pip:

```bash
$ pip install gcgc
```

There's also a docker image that can be downloaded that includes the necessary
dependencies to run Word Piece or BPE.

```bash
$ docker pull docker.io/thauck/gcgc
```

## Quick Start

```python
>>> from gcgc.tokenizer import kmer_tokenzier

>>> k_settings = kmer_tokenzier.KmerTokenizerSettings(alphabet='unambiguous_rna')
>>> k_tokenizer = kmer_tokenzier.KmerTokenizer(settings=k_settings)

>>> k_tokenizer.encode("AUGC")
[2, 3, 1, 4]
```

## General Approach

For the two tokenizers available in GCGC there is a settings class (inheriting
from `SequenceTokenizerSettings`) and a tokenizer class that preforms the work
(inheriting from `SequenceTokenizer`).

Altering settings can be done two ways. When the `SequenceTokenizer` derived
tokenizer is instantiated pass it a `SequenceTokenizerSettings` object. Or, set
the environment variables associated with the settings and they'll be used as
default.

Explicit creation:

```python
from gcgc.tokenizer import kmer_tokenzier

k_settings = kmer_tokenzier.KmerTokenizerSettings(alphabet="unambiguous_rna")
k_tokenizer = kmer_tokenzier.KmerTokenizer(settings=k_settings)
```

Via environment variables:

```python
import os
from gcgc.tokenizer import kmer_tokenzier

os.environ["GCGC_ALPHABET"] = "unambiguous_rna"

k_tokenizer = kmer_tokenzier.KmerTokenizer()
```

### Vocabulary

The vocabulary that maps token to integer used by GCGC is available at
`SequenceTokenizer.vocab`.

```python
# Using the same k_tokenizer from above.
>>> k_tokenizer.vocab
{'?': 0, 'G': 1, 'A': 2, 'U': 3, 'C': 4}
```

`?` is analogous to `UNK` in NLP tokenization, but chosen to maintain single
letter tokens.

### Finding the Settings for a Tokenizer

The settings objects are derived from `pydantic.BaseSettings` so calling
`.json_schema()` prints the schema. This also shows which environment variables
are associated with which settings.

```python
>>> print(kmer_tokenzier.KmerTokenizerSettings.schema_json(indent=True))
{ # only one option included for brevity
 "title": "KmerTokenizerSettings",
 "description": "The specification for the tokenizer.",
 "type": "object",
 "properties": {
  "bos_token": {
   "title": "Bos Token",
   "env": "GCGC_BOS_TOKEN",
   "env_names": [
    "gcgc_bos_token"
   ],
   "type": "string"
  }
 },
 "additionalProperties": false
}
```

## Available Tokenizers

GCGC implements two types of tokenizers:

- `KmerTokenizer` is a deterministic tokenizer that works similarly to n-grams.
- `BioSequencePiece` is a pre-training method that uses [SentencePiece](https://github.com/google/sentencepiece) by google.

## Links

- **Source Code** [GitHub Repo](https://github.com/tshauck/gcgc)

## Documentation Version

The documentation you're reading was build for version: `0.12.1`.
