GCGC doesn't doesn't directly integrate with any libraries, but it obviously
tries to make it easy for the developer to use their library of choice.

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

settings = gcgc.SequenceTokenizerSettings(alphabet="ATCG", kmer_size=2, kmer_step_size=2)
tokenizer = gcgc.SequenceTokenizer(settings=settings)

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

