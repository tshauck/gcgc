# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""CLI Functions for GCGC."""

import json
import pathlib

import click
from Bio import SeqIO

from gcgc.tokenizer.sentence_piece_tokenizer import BioSequencePiece, BioSequencePieceSettings
from gcgc.tokenizer.kmer_tokenzier import KmerTokenizer, KmerTokenizerSettings


@click.group()
def cli():
    """Console script for gcgc."""


@cli.command()
def version():
    """Print the version and exit."""
    click.echo("0.12.0-dev.5")


@cli.group("tokenizer")
def tokenizer_group():
    """Entrypoint for the tokenizer command."""


@tokenizer_group.command("kmer")
@click.argument("input_fasta", type=str)
@click.option(
    "--vocab_path", type=str, default="kmer_vocab.json", help="Where to save the vocabulary."
)
def kmer_tokernizer(input_fasta: str, vocab_path: str):
    """Run the kmer tokenizer.

    \b
    INPUT_FASTA is the path to the FASTA file to use for training.

    """
    tokenizer = KmerTokenizer(settings=KmerTokenizerSettings())
    with open(input_fasta) as input_handler:
        for record in SeqIO.parse(input_handler, format="fasta"):
            tokenized_sequence = tokenizer.encode(str(record.seq))
            output_record = {
                "token_ids": tokenized_sequence,
                "sequence_id": record.id,
                "n_tokens": len(tokenized_sequence),
            }
            click.echo(json.dumps(output_record))

    if vocab_path:
        with open(vocab_path, "w") as vocab_file:
            json.dump(tokenizer.vocab, vocab_file, indent=2)


@tokenizer_group.command("sentencepiece")
@click.argument("input_fasta", type=str)
@click.argument("model_prefix", type=str)
def sentencepiece_tokenizer(input_fasta: str, model_prefix: str):
    """Use a pretrained sentencepiece tokenizer to separate a FASTA into a space separate text file.

    See `gcgc train sentencepiece` for information on how to generate the pretrained model.

    \b
    INPUT_FASTA is the path to the FASTA file to use for training.
    MODEL_PREFIX is the path of the pretrained model.

    """
    tokenizer = BioSequencePiece(settings=BioSequencePieceSettings(model_prefix=model_prefix))

    with open(input_fasta) as input_handler:
        for record in SeqIO.parse(input_handler, format="fasta"):
            tokenized_sequence = tokenizer.encode(str(record.seq))
            output_record = {
                "token_ids": tokenized_sequence,
                "sequence_id": record.id,
                "n_tokens": len(tokenized_sequence),
            }
            click.echo(json.dumps(output_record))


@cli.group()
def train():
    """Entrypoint for the train command."""


@train.command("sentencepiece")
@click.argument("input_fasta", type=str)
@click.argument("model_prefix", type=str)
def sentencepiece_train(input_fasta: str, model_prefix: str):
    """Use sentencepiece to fit a subword tokenization model.

    See "SentencePiece: A simple and language independent subword tokenizer and detokenizer for
    Neural Text Processing" (Kudo 2018) for more information
    (https://github.com/google/sentencepiece).

    \b
    INPUT_FASTA is the path to the FASTA file to use for training.
    MODEL_PREFIX is the prefix to use for the trained model and the vocabulary.

    """
    input_path = pathlib.Path(input_fasta)
    if not input_path.exists():
        raise ValueError(f"{input_path} does not exist.")

    bsp = BioSequencePiece(settings=BioSequencePieceSettings(model_prefix=model_prefix))
    bsp.fit_on_fasta(input_path)
