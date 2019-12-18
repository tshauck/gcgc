# GCGC

GCGC is a package written in Python for pre-processing for biological sequences. Think of it like a
Natural Language Processing pre-processing toolkit with design choices oriented towards the
differences the sequences found in natural language vs biology.

At its core is a tokenizer specifically for biological sequences that can be
used as part of PyTorch or TensorFlow workflows.

<!-- vim-markdown-toc GFM -->

* [Installation](#installation)
* [Using the Tokenizer](#using-the-tokenizer)
* [Using with `torchtext`](#using-with-torchtext)
* [Links](#links)
* [Documentation Version](#documentation-version)

<!-- vim-markdown-toc -->

## Installation

Install GCGC via pip:

```bash
$ pip install gcgc
```

## Using the Tokenizer

To use the tokenizer first instantiate the object, then that object can be
called to perform the tokenization.

```python
import gcgc

spec = gcgc.SequenceTokenizerSpec(alphabet="ATCG")
tokenizer = gcgc.SequenceTokenizer(spec)

sequence = "ATCG"
output_sequence = tokenizer(sequence)

assert output_sequence == ["A", "T", "C", "G"]
```

Like you see, by default, the tokenizer will tokenize a kmer of size 1 and take
a step of size 1, essentially converting the string to a list.

By altering the `kmer_size` and `kmer_step_size` parameters passed to
`SequenceTokenizerSpec`, how the sequence is broken into kmers can be controlled.

`kmer_size` control the size of the kmer, so if the desire is to break the
sequence into 2mers, then pass `kmer_size=2`.

```python
import gcgc

spec = gcgc.SequenceTokenizerSpec(alphabet="ATCG", kmer_size=2)
tokenizer = gcgc.SequenceTokenizer(spec)

sequence = "ATCG"
output_sequence = tokenizer(sequence)

assert output_sequence == ["AT", "TC", "CG"]
```

Then `kmer_step_size` controls how many kmers to pass through each step. In the
example with `kmer_size=2` the output produced three tokens because while the
window size (`kmer_size`) was two, there was an output token at each step. If
`kmer_step_size=2` too, then there will be two output tokens.

```python
from gcgc import SequenceTokenizerSpec

spec = SequenceTokenizerSpec(alphabet="ATCG", kmer_size=2, kmer_step_size=2)
tokenizer = gcgc.SequenceTokenizer(spec)

sequence = "ATCG"
output_sequence = tokenizer(sequence)

assert output_sequence == ["AT", "CG"]
```

## Using with `torchtext`

The tokenizer implemented in `gcgc` is compatible with `torchtext`'s `Field`
class.

So imagine there's a TSV file that has a sequence of interest. Generally there'd
also be other metadata, like a label.

```
sequence
ATCG
ATTT
ACGG
```

It can be loaded using `torchtext`'s `TabularDataset` object, but instead of
supplying no tokenize function to `Field`, we can supply GCGC's tokenizer.

```python
import gcgc
from torchtext import data

spec = gcgc.SequenceTokenizerSpec(alphabet="ATCG", kmer_size=2, kmer_step_size=2)
tokenizer = gcgc.SequenceTokenizer(spec)

sequence_field = data.Field(
    tokenize=tokenizer,
    batch_first=True,
    fix_length=10
)

train = data.TabularDataset(
    path="./data.tsv",
    format="tsv",
    skip_header=True,
    fields=[("sequence", sequence_field)],
)

sequence_field.build_vocab(train, min_freq=1)

train_iter = data.Iterator(
    train,
    sort_key=lambda x: len(x.sequence),
    batch_size=32,
    sort_within_batch=True,
    repeat=False,
)

for example in train_iter:
    print(example.sequence)
    break
```

We'd then get a tensor for the sequence attribute on the example:

```python
tensor([[2, 4, 1, 1, 1, 1, 1, 1, 1, 1],
        [3, 5, 1, 1, 1, 1, 1, 1, 1, 1],
        [2, 6, 1, 1, 1, 1, 1, 1, 1, 1]])
```

## Links

- **Source Code** [GitHub Repo](https://github.com/tshauck/gcgc)

## Documentation Version

The documentation you're reading was build for version: `0.12.0-dev.2`.
