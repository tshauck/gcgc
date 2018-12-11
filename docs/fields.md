# Fields

Fields are additional pieces of data, including labels (the thing being predicted).

At a most basic level, a field has a name, and depending on the type of field it may take additional
options to do more work.

Once a field is declared the it is added to the data dictionary that is returned when processing
sequences. See the example in the PyTorch integration section.

## Types of Fields

As stated different kinds of fields afford different functionality.

### `LabelField`

A `LabelField` ascribes a particular label to a sequence.

### `FileMetaDataField`

A `FileMetaDataField` can be used in cases where the filenames contain label information.

For example, imagine there's a list of files:

```python
from pathlib import Path
paths = [Path("a/f1.fasta"), Path("a/f2.fasta"), Path("b/f3.fasta")]
```

In this case there are two files, f1 and f2, that are under label "a" and f3 is under label b. In
order to create a `FileMetaDataField` we'll write a function that takes in a `Path` and returns a
string.

```python
def preprocess(p: Path) -> str:
    # Path("a/b") -> Path("a") -> "a"
    return str(p.parent)
```

Then create the `FileMetaDataField`.

```python
from gcgc.fields import FileMetaDataField
label = FileMetaDataField.from_paths(name='label', paths=paths, preprocess=preprocess)
```

After that, the field has an encoding dictionary and, more importantly, it will add the label name and
corresponding integer encoded class to the dataset.

```python
label.encoding_dict
{'a': 0, 'b': 1}
```
