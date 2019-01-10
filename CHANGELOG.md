# Change Log

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
