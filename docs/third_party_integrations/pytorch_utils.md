# PyTorch Integration

GCGC implements a `torch.utils.data.Dataset` that loads data from common bioinformatics file
formats via BioPython.

The easiest way to see how this works to see [Splice Site](./../examples/splice-site.md) example.

## Genomic Dataset

The `GenomicDataset` implements the interface for a PyTorch `Dataset`.

Assuming you have the FASTAs to read, create a dataset by first creating a
`TorchSequenceParser` object, the alphabet, then the actual `Dataset`.

```python
from gcgc.alphabet import IUPACProteinEncoding
from gcgc.ml.pytorch_utils.data import GenomicDataset
from gcgc.ml.pytorch_utils.parser import TorchSequenceParser

length_parser = EncodedSeqLengthParser(conform_to=100)  # Make seq length uniform at 100
parser = TorchSequenceParser(
  encapsulate=False, seq_length_parser=length_parser, sequence_offset=1
)

alphabet = IUPACProteinEncoding()  # Use an Amino Acid vocab.

files = ['fasta1.fasta', 'fasta2.fasta']
dataset = GenomicDataset.from_paths(files, parser, alphabet=alphabet)
```

- `GenomicDataset`, again, implements the interface.
- `TorchSequenceParser` is an object that helps parse the incoming sequences to
  PyTorch tensors.
- `IUPACProteinEncoding` is the protein encoding alphabet, think of it as the
  vocabulary from NLP.

Now that a `Dataset` exists, it can be passed to `DataLoader`.

```python
from torch.utils.data import DataLoader
data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

for batch in data_loader:

  # seq_tensor is (B, S) tensor where B is the batch size and A is the sequence
  # length, which in this case is 100
  seq_tensor = batch['seq_tensor']
```

### Use with PyTorch LSTM Models

One important gotcha is that for using PyTorch models is that typically it
expects the batches to have the sequence length as the first dimension and batch
size as the second. Therefore, you may want to transpose the matrix.

```python
# Now has batch first.
seq_tensor = seq_tensor.transpose(0, 1)
```
