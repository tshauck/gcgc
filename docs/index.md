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
from gcgc.tokenizer import SequenceTokenizer

tokenizer = SequenceTokenizer()

sequence = "ATCG"
output_sequence = tokenizer(sequence)

assert output_sequence == ["A", "T", "C", "G"]
```

Like you see, by default, the tokenizer will tokenize a kmer of size 1 and take
a step of size 1, essentially converting the string to a list.

By altering the `kmer_size` and `kmer_step_size` parameters passed to
`SequenceTokenizer`, how the sequence is broken into kmers can be controlled.

`kmer_size` control the size of the kmer, so if the desire is to break the
sequence into 2mers, then pass `kmer_size=2`.

```python
from gcgc.tokenizer import SequenceTokenizer

tokenizer = SequenceTokenizer(kmer_size=2)

sequence = "ATCG"
output_sequence = tokenizer(sequence)

assert output_sequence == ["AT", "TC", "CG"]
```

Then `kmer_step_size` controls how many kmers to pass through each step. In the
example with `kmer_size=2` the output produced three tokens because while the
window size (`kmer_size`) was two, there was an output token at each step. If
`kmer_step_size=2` too, then there will be two output tokens.

```python
from gcgc import SequenceTokenizer

tokenizer = SequenceTokenizer(kmer_size=2, kmer_step_size=2)

sequence = "ATCG"
output_sequence = tokenizer(sequence)

assert output_sequence == ["AT", "CG"]
```

## Using with `torchtext`

The tokenizer implemented in `gcgc` is compatible with `torchtext`'s `Field`
class.

So imagine there's a TSV file that has a label of interest and a sequence.

```
label   sequence
A   ATCG
B   ATTT
A   ACGG
```

It can be loaded using `torchtext`'s `TabularDataset` object, but instead of
supplying no tokenize function to `Field`, we can supply GCGC's tokenizer.

```python
from gcgc import SequenceTokenizer
from torchtext import data

dna_tokenizer = SequenceTokenizer(kmer_size=2, kmer_step_size=2)

sequence_field = data.Field(
    tokenize=dna_tokenizer,
    batch_first=True,
    fix_length=40,
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
tensor([[13,  2, 13,  ...,  1,  1,  1],
        [12,  8,  3,  ...,  1,  1,  1],
        [12,  7,  5,  ...,  1,  1,  1],
        ...,
        [ 8,  8, 12,  ...,  1,  1,  1],
        [17, 14,  9,  ...,  1,  1,  1],
        [ 7,  8,  9,  ...,  1,  1,  1]])
```

The `sequence_field` will now have the relevant vocabulary having used GCGC's
tokenizer.

```
import pprint
print(pprint.pformat(sequence_field.vocab.stoi))
defaultdict(<bound method Vocab._default_unk_index of <torchtext.vocab.Vocab object at 0x1243bca90>>,
            {'<pad>': 1,
             '<unk>': 0,
             'AA': 12,
             'AC': 14,
             'AG': 4,
             'AN': 21,
             'AS': 22,
             'AT': 15,
             'CA': 7,
             'CC': 2,
             'CG': 17,
             'CT': 3,
             'DG': 23,
             'GA': 8,
             'GC': 9,
             'GD': 24,
             'GG': 5,
             'GN': 25,
             'GT': 13,
             'NA': 26,
             'NC': 27,
             'NG': 28,
             'NN': 18,
             'NT': 20,
             'RT': 29,
             'TA': 16,
             'TC': 11,
             'TG': 6,
             'TN': 19,
             'TT': 10}
```

## Links

- **Source Code** [GitHub Repo](https://github.com/tshauck/gcgc)

## Documentation Version

The documentation you're reading was build for version: `0.11.0-dev.3`.
