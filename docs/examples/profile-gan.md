# Profile GAN

In this example, we build a proof-of-concept Profile GAN drawing inspiration from the idea of a
Profile HMM.

## Profile HMM

A Profile HMM is a Hidden Markov Model that allows the user to fit a model from a profile of
sequences for the purposes of alignment or other sequence search tasks.

There is a wealth of information available on [Profile HMMs](https://en.wikipedia.org/wiki/HMMER#Profile_HMMs).

## Profile GAN

With a HMM the sequence likelihood is tested against a random model to tell if it's likely the model
generated the sequence or a random model. In this case, we'll use the adverserial aspects of the GAN
to train a generator and a discriminator, with the goal the generator can generate sequences like
the profile sequences and the discriminator can be used to score sequences as being part of the
profile or not.

## Acquiring Data

As an example, we'll use the Bacterial-like globin family, which can be downloaded from [EMBL-EBI](https://pfam.xfam.org/family/PF01152#tabview=tab0).

## Loading the Data

```python
from gcgc.ml.pytorch_utils.parser import TorchSequenceParser

parser = TorchSequenceParser()

from gcgc.alphabet import ExtendedIUPACProteinEncoding
alphabet = ExtendedIUPACProteinEncoding(add_lower_case_for_inserts=True)

from gcgc.ml.pytorch_utils.data import GenomicDataset
from gcgc.tests.fixtures import PF01152_PATH_FULL

dataset = GenomicDataset.from_path(PF01152_PATH_FULL, parser, alphabet=alphabet)

from torch.utils.data import DataLoader

data_loader = DataLoader(dataset, batch_size=128, shuffle=True)

for i in range(10):  # 10 epochs
    for batch in data_loader:  # Created from the GenomicDataset
        sequences = batch["seq_tensor"].to("cpu")
```
