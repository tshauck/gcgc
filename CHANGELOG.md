# Changelog

## 0.12.0 (Unreleased)

- Improved the docs to reflect the `SequenceTokenizerSpec` that was added in
  0.11.0.
- Made max length optional for the tokenizer.
- Added CLI that parses use the Sequence Piece library.
- Began versioning docker build, and make pushing easier.
- Have the tokenizer resolve the named alphabets.

## 0.11.0 (2019-11-15)

### Added

- Added the `SequenceTokenizerSpec` object for specifying the tokenizer.
- Added `Vocab` object for storing the int to token, and token to int encodings.
- Added example of using tensorflow/keras together with gcgc.

## 0.10.0 (2019-11-09)

### Changed

`gcgc` has been revamped quite a bit to better support existing processing
pipelines for NLP without trying to do to much. See the docs for more
information about how this works.

## 0.9.0 (2019-08-05)

### Added

- Parser now outputs the length of the tensor not including padding. This is
  useful for packing and length based iteration.
- Generating masked output from the parse_record method is now available.
- Alphabet can include an optional mask token.

### Changed

- Can now specify how large of kmer step size to generate when supplying a kmer
  value.
- Renames EncodedSeq.integer_encoded to EncodedSeq.get_integer_encoding which
  takes a kmer_step_size to specify how large of steps to take when encoding.
- Add parsed_seq_len to the SequenceParser object to control how much padding to
  apply to the end of the integer encoded sequence. This is useful since a batch
  of tensors is expected to have the same size.

## 0.8.0 (2019-07-04)

### Fixed

- Broken test due to platform differences in `Path.glob` sorting.

### Added

- User can specify to use start or end tokens optionally.

### Removed

- Removed one_hot_encoding. The user can do that pretty easily if needed. E.g.
  see `scatter` in PyTorch.

## 0.7.0 (2019-06-22)

### Added

- Properties to access the integer encodings of special tokens. (35cae2a)
  - `Alphabet.encoded_start`
  - `Alphabet.encoded_end`
  - `Alphabet.encoded_padding`
- Remove uniprot dataset creation. (e233162)
- Simplify index handling for GenomicDataset. (3213a9e)

## 0.6.1 (2019-06-10)

### Added

- Updated package management so gcgc is easier to use with other version of
  torch.

## 0.6.0 (2019-04-04)

### Added

- Ability for kmer size to be passed to an alphabet.

## 0.5.2 (2019-03-21)

### Added

- Add Dockerfile and docker-compose.yml for development.
- `EncodedSeq.shift`, which will shift sequence by an offset integer.
- `EncodedSeq.from_integer_encoded_seq` will take a list of integers and an
  alphabet and return an EncodedSeq object.
- Add the ability to apply a function to the rollout_kmers yielded values.

### Changed

- Alphabet special characters are now located at the start, rather than the end,
  of the letters and token sequence.

## 0.5.1 (2019-01-09)

### Added

- Add extra css to make underline links in articles.
- Exit if the download directory doesn't exist in the call to download organism.
- Wording improvements in docs.

## 0.5.0 (2018-12-31)

### Added

- Include `seq_tensor_one_hot` in the PyTorch Parser.
- Added a `GCGCRecord.encoded_seq` property.
- New `gcgc.random` module to start holding sequence data.
- New `gcgc.rollout` module to handle working through chunks of sequences.
  - `rollout_kmers` will roll out [kmers][1].
  - `rollout_seq_features` will roll out the `SeqFeatures` from a `SeqRecord`.
- `EncodingAlphabet` now can optionally take a `gap_characters` set of characters to add to the
  alphabet letters. It also takes `add_lower_case_for_inserts` which will duplicate the alphabet,
  but convert the letters to lowercase.

### Changed

### Fixed

- Fixed bug in `GenomicDataset.from_path` where it still referred to `init_from_path_generator`.

## 0.4.0

### Added

- `EncodedSeq` now supports iterating through kmers, see `EncodedSeq.rollout_kmers` for options.
- GCGC is citable.
- GCGC now has a CHANGELOG.md.

[1]: https://en.wikipedia.org/wiki/K-mer
