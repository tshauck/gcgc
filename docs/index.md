# Welcome to the GCGC Docs

GCGC is a package written in Python to facilitate feature pre-processing for biological sequences.
Think of it like a Natural Language Processing pre-processing toolkit with design choices oriented
towards the differences between NLP and bio sequences.

There are two main points of entry to GCGC. First, it can be imported as a Python package, and,
second, used as a command line tool for specific types of transformations. In this case, the ideas
is for the tools to be integrated into a larger processing pipeline.

## Python Package

In general, operations that would modify the sequence return a new `EncodedSeq` object, which is a subclass of `Bio.Seq.Seq`.

```python
from gcgc.encoded_seq import EncodedSeq
from gcgc.alphabet import UnambiguousDnaAlphabet

>>> es = EncodedSeq("ATCG", UnambiguousDnaAlphabet())
EncodedSeq('ATCG', UnambiguousDnaAlphabet())

>>> es.pad(pad_to=10)
EncodedSeq('ATCG||||||', UnambiguousDnaAlphabet())

>>> es.encapsulate()
EncodedSeq('>ATCG<', UnambiguousDnaAlphabet())
```

Otherwise there exists properties on the `EncodedSeq` object to get it in a form amenable to traditional ML modeling pipelines.

For example, turning it into a one hot encoded matrix.

```python
>>> es.one_hot_encode_sequence
[[0, 1, 0, 0, 0, 0, 0],
 [0, 0, 1, 0, 0, 0, 0],
 [0, 0, 0, 1, 0, 0, 0],
 [1, 0, 0, 0, 0, 0, 0]]
```

To do this, the alphabet holds a mapping from the letters on the associated `Bio.Alphabet.Alphabet` object.

```python
>>> es.alphabet.encoding_index
{'G': 0, 'A': 1, 'T': 2, 'C': 3, '>': 4, '<': 5, '|': 6}
```

The last three characters are the start (`es.alphabet.START`), the end (`es.alphabet.END`), and the padding (`es.alphabet.PADDING`), respectively.

The alphabet concept is once of the key differences relative to NLP vocabulary, the in the bio case the vocabulary is relatively small and known a priori.

## Command Line Interface
