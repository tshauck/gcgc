# Sequence Parsers

In GCGC, parsers specify how the input data goes from the raw BioPython SeqRecords into data
suitable for ML. They act as a higher level abstraction above Alphabets and EncodedSeqs.

## Third Party Parsers

SequenceParser's subclasses also act as the intermediary between the simple data types in GCGC and
the specific type native to the library in question.

For example, to generate PyTorch tensors, use the `TorchSequenceParser`.
