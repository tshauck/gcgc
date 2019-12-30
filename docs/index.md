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

* `KmerTokenizer` is a deterministic tokenizer that works similarly to n-grams.
* `BioSequencePiece` is a pre-training method that uses [SentencePiece](https://github.com/google/sentencepiece) by google.

## Using the CLI

As mentioned, GCGC can be used from via a container, but regardless, when
installed with pip the `gcgc` command is available.

To get help:

```
gcgc -h
```

Though generally `gcgc` expects to operate on a FASTA file with options passed via env vars.

```console
$ GCGC_ALPHABET=protein gcgc tokenizer kmer sprot_20k.fasta | head -n 1
{"token_ids": [11, 1, 5, 16, 1, 4, 3, 18, 10, 9, 4, 20, 3, 15, 15, 15, 15, 11, 4, 1, 10, 10, 10, 16, 10, 20, 20, 13, 12, 3, 15, 9, 10, 10, 3, 20, 9, 4, 19, 16, 13, 13, 15, 18, 14, 18, 4, 2, 13, 9, 1, 13, 18, 4, 19, 12, 12, 13, 13, 16, 4, 9, 6, 10, 8, 18, 6, 7, 5, 16, 6, 8, 9, 20, 9, 6, 4, 9, 1, 14, 1, 16, 4, 18, 3, 18, 12, 9, 11, 2, 2, 19, 18, 16, 9, 5, 9, 3, 1, 11, 15, 15, 20, 14, 6, 8, 14, 17, 2, 9, 8, 13, 6, 9, 18, 10, 16, 3, 10, 3, 1, 9, 8, 9, 1, 20, 12, 10, 17, 18, 4, 6, 18, 4, 6, 5, 18, 15, 20, 16, 15, 18, 17, 9, 14, 7, 18, 1, 1, 5, 10, 9, 4, 10, 15, 7, 16, 9, 14, 20, 4, 12, 18, 12, 10, 8, 7, 20, 8, 10, 17, 3, 9, 15, 18, 3, 8, 14, 7, 10, 4, 9, 3, 10, 18, 9, 3, 5, 9, 1, 10, 18, 4, 16, 1, 7, 15, 11, 15, 14, 6, 7, 11, 8, 12, 18, 9, 20, 8, 10, 20, 14, 10, 10, 9, 9, 7, 6, 7, 6, 13, 3, 6, 13, 3, 8, 10, 17, 18, 9, 17, 6, 16, 9, 6, 18, 10, 20, 3, 3, 16, 5, 15, 9, 8, 20, 17, 3, 10, 6, 19, 9, 5, 17, 13, 10], "sequence_id": "sp|Q6GZX4|001R_FRG3G", "sequence_description": "sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1", "n_tokens": 256}
{"token_ids": [11, 16, 8, 8, 6, 1, 17, 15, 10, 14, 12, 3, 9, 16, 3, 17, 20, 16, 1, 6, 13, 2, 20, 1, 6, 6, 2, 16, 1, 5, 17, 13, 15, 6, 17, 2, 6, 9, 3, 19, 3, 10, 6, 4, 14, 17, 2, 1, 16, 6, 5, 2, 17, 16, 14, 13, 10, 2, 1, 15, 8, 9, 9, 17, 14, 18, 2, 6, 10, 15, 20, 16, 16, 9, 6, 9, 3, 13, 10, 18, 16, 1, 4, 19, 3, 16, 15, 6, 1, 13, 20, 18, 15, 2, 17, 20, 3, 1, 3, 10, 8, 3, 17, 14, 1, 14, 18, 3, 14, 5, 18, 16, 11, 5, 6, 4, 16, 13, 16, 10, 1, 4, 15, 20, 2, 11, 15, 6, 18, 9, 12, 17, 1, 6, 4, 10, 18, 16, 15, 18, 16, 16, 3, 1, 3, 13, 1, 6, 6, 19, 2, 15, 9, 19, 20, 16, 1, 7, 15, 6, 13, 3, 14, 3, 1, 1, 10, 6, 16, 5, 2, 8, 9, 12, 13, 6, 1, 1, 3, 2, 9, 2, 8, 12, 15, 1, 16, 3, 13, 18, 20, 14, 9, 18, 9, 17, 10, 7, 1, 20, 13, 3, 14, 2, 19, 20, 18, 13, 2, 1, 1, 3, 18, 6, 4, 10, 9, 11, 6, 17, 14, 15, 3, 17, 13, 17, 12, 2, 13, 17, 14, 18, 2, 14, 8, 18, 5, 12, 11, 10, 3, 3, 6, 16, 18, 17, 11, 3, 3, 18, 9, 12, 17, 8, 12, 2, 3, 5, 16, 9, 20, 18, 13, 13, 13, 13, 13, 13, 9, 13, 17, 13, 13, 17, 13, 13, 17, 13, 13, 17, 13, 13, 17, 13, 13, 17, 13, 13, 17, 13, 13, 17, 13, 15, 13, 18, 7, 12, 15, 9, 18, 11, 5, 5, 18, 1, 6, 1, 18, 10, 18, 1, 8, 10, 8, 16, 17, 18, 15, 19], "sequence_id": "sp|Q6GZX3|002L_FRG3G", "sequence_description": "sp|Q6GZX3|002L_FRG3G Uncharacterized protein 002L OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-002L PE=4 SV=1", "n_tokens": 320}
```

### Training the `BioSequencePiece` Tokenizer

The `BioSequencePiece` can be trained via the command line:

```console
$ gcgc train sentencepiece ./gcgc/data/sprot_20k.fasta swissprot_tokenizer
sentencepiece_trainer.cc(116) LOG(INFO) Running command: --input=/var/folders/93/cn514jtx3_5djq73slspkjr40000gn/T/tmpxpu5mmjb/input_textfiles.txt --model_prefix=swissprot_tokenizer --vo
cab_size=8000 --model_type=unigram --max_sentence_length=4192 --unk_piece=? --bos_id=-1 --eos_id=-1 --pad_id=-1
sentencepiece_trainer.cc(49) LOG(INFO) Starts training with :
TrainerSpec {
  input: /var/folders/93/cn514jtx3_5djq73slspkjr40000gn/T/tmpxpu5mmjb/input_textfiles.txt
  input_format:
  model_prefix: swissprot_tokenizer
  model_type: UNIGRAM
  vocab_size: 8000
  self_test_sample_size: 0
  character_coverage: 0.9995
  input_sentence_size: 0
  shuffle_input_sentence: 1
  seed_sentencepiece_size: 1000000
  shrinking_factor: 0.75
  max_sentence_length: 4192
  num_threads: 16
  num_sub_iterations: 2
  max_sentencepiece_length: 16
  split_by_unicode_script: 1
  split_by_number: 1
  split_by_whitespace: 1
  treat_whitespace_as_suffix: 0
  hard_vocab_limit: 1
  use_all_vocab: 0
  unk_id: 0
  bos_id: -1
  eos_id: -1
  pad_id: -1
  unk_piece: ?
  bos_piece: <s>
  eos_piece: </s>
  pad_piece: <pad>
  unk_surface:  ⁇
}
...
```

After this runs, the model and vocabulary are available in the working
directory.

```console
$ ls swissprot_tokenizer.*
swissprot_tokenizer.model  swissprot_tokenizer.vocab
```

## Links

- **Source Code** [GitHub Repo](https://github.com/tshauck/gcgc)

## Documentation Version

The documentation you're reading was build for version: `0.12.0-dev.6`.
