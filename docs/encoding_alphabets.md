# Encoding Alphabets

The concept of an Alphabet is an analog to the concept of a vocabulary in NLP, i.e. the set of
possible tokens present in the corpus, or in the biological case, the sequences.

For example, in the case of DNA where there is no ambiguity in the reads, the alphabet is `{A, T, C, G}`.

## Implementation

Like other areas of GCGC, BioPython is used to supply the basics of the Alphabet with GCGC adding
additional functionality of the ML case.

In general, for a BioPython Alphabet, there exists a corresponding GCGC alphabet that handles
encoding and decoding a sequence into an integer or one-hot encoding.

For example, in the case described in the into the `IUPACUnambiguousDNAEncoding` alphabet would be a
reasonable choice.

```python
from gcgc.alphabet import IUPACUnambiguousDNAEncoding

alpha = IUPACAmbiguousDNAEncoding()
```

Once the alphabet is instantiated, it's simple to encode or decode values to go between text and
integer encodings.

```python
alpha.encode_token("A")
# 1

alpha.decode_token(1)
# "A"
```

## Custom Alphabets

Currently GCGC does not provide direct support for user defined Alphabets, though it is on the
roadmap. Hopefully, one of the standard [IUPAC](https://iupac.org/)'s Alphabets will work. See
the module `gcgc.alphabet.iupac` for more info.
