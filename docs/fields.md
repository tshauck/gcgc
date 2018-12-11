# Fields

Fields are additional pieces of data, including labels (the thing being predicted).

At a most basic level, a field has a name, and depending on the type of field it may take additional
options to do more work.

## Types of Fields

As stated different kinds of fields afford different functionality.

### `LabelField`

A `LabelField` is used to ascribe a particular label to a sequence.

### `FileMetaDataField`

A `FileMetaDataField` can be used in cases where the file names contain labels.
