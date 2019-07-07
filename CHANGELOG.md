# Changelog

## Development

### Added

- User can specify to use start or end tokens optionally.
- Fixed bug in one_hot_encoded property due to dimensionality difference

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
