# CLI

GCGC (as of version 0.12) includes a command line interface to facilitate
tokenization file directly.

## Basic Use

Assuming `gcgc` is installed, the best first thing to do is to look at the help.

```console
$ gcgc --help
Usage: gcgc [OPTIONS] COMMAND [ARGS]...

  Console script for gcgc.

Options:
  --help  Show this message and exit.

Commands:
  tokenizer  Entrypoint for the tokenizer command.
  train      Entrypoint for the train command.
  version    Print the version and exit.
```

The two important commands are:

- tokenizer: which will perform tokenization of an input file.
- train: which is used for training tokenizers, assume they're needed

### Tokenizer

The tokenizer will take an input file (FASTA is currently the only file type
supported) and convert it into the JSON newline delimited format.

As a realistic example, this will use the `protein` alphabet to tokenize the FASTA
found at `./uniprot-reviewed_yes+taxonomy_83333.fasta` and stream the data into
`swissprot-83333.jsonl`. It also outputs the vocabulary that the tokenizer used
to `kmer-vocab.json`.

```bash
$ GCGC_ALPHABET=protein gcgc tokenizer kmer \
  --vocab_path ./kmer-vocab.json \
  ./uniprot-reviewed_yes+taxonomy_83333.fasta > swissprot-83333.jsonl
```

Looking at the first row we get:

```console
$ head -n 1 ./swissprot-83333.jsonl | jq
{
  "token_ids": [
    11,
    12,
    12  # omitted the rest of the tokens.
  ]
  "sequence_id": "sp|P0AFL1|POTI_ECOLI",
  "sequence_description": "sp|P0AFL1|POTI_ECOLI Putrescine transport system permease protein PotI OS=Escherichia coli (strain K12) OX=83333 GN=potI PE=3 SV=1",
  "n_tokens": 281
}
```

And the vocab:

```console
cat ./kmer-vocab.json
{
  "?": 0,
  "A": 1,
  "C": 2,
  "D": 3,
  "E": 4,
  "F": 5,
  "G": 6,
  "H": 7,
  "I": 8,
  "K": 9,
  "L": 10,
  "M": 11,
  "N": 12,
  "P": 13,
  "Q": 14,
  "R": 15,
  "S": 16,
  "T": 17,
  "V": 18,
  "W": 19,
  "Y": 20
}
```

Generally options to the CLI are limited to those supporting I/O, and internal
options such as alphabet or stride length are expected to be set using
environment variables.

### Train

The train command is only relevant to tokenizers that require use the underlying
data for parameter estimation before the actual tokenization process. For GCGC,
that means sentencepiece.

Using the same dataset, the tokenizer is trained first, then used second.

```bash
$ gcgc train sentencepiece ./uniprot-reviewed_yes+taxonomy_83333.fasta ./swissprot-sp-83333
```

This will use GCGC's bindings to Google's sentencepiece package to facilitate
training on the underlying FASTA.

The output of vocab and serialized model is stored at the prefix `./swissprot-sp-83333.{vocab,model}`.

Then once the tokenizer is trained, it can be applied to the dataset producing
the same output structure as the deterministic tokenizer.

```console
$ gcgc tokenizer sentencepiece ./uniprot-reviewed_yes+taxonomy_83333.fasta ./swissprot-sp-83333 | head -n 1 | jq
{
  "token_ids": [
    2883,
    526,
    51,
    2915,
    7714,
    2975,
    2212,
    2924,
    30,
    2003,
    3014
  ]
  "sequence_id": "sp|P0AFL1|POTI_ECOLI",
  "sequence_description": "sp|P0AFL1|POTI_ECOLI Putrescine transport system permease protein PotI OS=Escherichia coli (strain K12) OX=83333 GN=potI PE=3 SV=1",
  "n_tokens": 97
}
```

## Docker

To help facilitate using GCGC as part of larger bioinformatic or DL workflows,
gcgc publish an image to docker.io which can be pull like other images. The one
difference is there is no latest, so you'll need to specify the version as the
tag.

```console
$ docker pull docker.io/thauck/gcgc:0.12.0
0.12.0: Pulling from thauck/gcgc
Digest sha256:f94f966c00d065a7c9c9ce52ec10c227c13b34211d3de76f08616162d63301ef
Status: Downloaded newer image for thauck/gcgc:0.12.0
docker.io/thauck/gcgc:0.12.0
```

Once the image is pulled, use the docker image like you would the regular CLI.

```bash
$ docker run --rm \
  -e GCGC_SP_MODEL_TYPE=unigram \
  -v $PWD:/workspace \
  -v $PWD/uniprot-reviewed_yes+taxonomy_83333.fasta:/data/uniprot-reviewed_yes+taxonomy_83333.fasta \
  -t docker.io/thauck/gcgc:0.12.0 \
  train sentencepiece /data/uniprot-reviewed_yes+taxonomy_83333.fasta /workspace/docker-sp
```

The complicating factor is the volume mounts. GCGC works in the `/workspace`
directory in the container, so it should be mounted such so it the host machine
-- here the host current working directory.

Given the trained model, the next step is to use the model for tokenization.

```console
$ docker run --rm \
  -v $PWD:/workspace \
  -v $PWD/uniprot-reviewed_yes+taxonomy_83333.fasta:/data/uniprot-reviewed_yes+taxonomy_83333.fasta \
  -t docker.io/thauck/gcgc:0.12.0\
  tokenizer sentencepiece /data/uniprot-reviewed_yes+taxonomy_83333.fasta /workspace/docker-sp \
  | head -n 1
{"token_ids": [2880, 533, 44, 2904, 7702], "sequence_id": "sp|P0AFL1|POTI_ECOLI", "sequence_description": "sp|P0AFL1|POTI_ECOLI Putrescine transport system permease protein PotI OS=Escherichia coli (strain K12) OX=83333 GN=potI PE=3 SV=1", "n_tokens": 97}
```
