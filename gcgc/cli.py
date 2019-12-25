# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""CLI Functions for GCGC."""

import json
import logging
import pathlib
import sys

import click
from Bio import SeqIO

import sentencepiece as spm

from gcgc import __version__
from gcgc.tokenizer.sentence_piece_tokenizer import BioSequencePiece, BioSequencePieceSettings

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
LOG = logging.getLogger(__name__)


@click.group()
def cli():
    """Console script for gcgc."""


@cli.command()
def version():
    """Print the version and exit."""
    print(__version__)


@cli.group()
def tokenizer():
    """Entrypoint for the tokenizer command."""


@tokenizer.command("sentencepiece")
@click.argument("input_fasta", type=str)
@click.argument("model_path", type=str)
def sentencepiece_tokenizer(input_fasta: str, model_path: str):
    """Use a pretrained sentencepiece tokenizer to separate a FASTA into a space separate text file.

    See `gcgc train sentencepiece` for information on how to generate the pretrained model.

    \b
    INPUT_FASTA is the path to the FASTA file to use for training.
    MODEL_PATH is the path of the pretrained model.

    """
    sentence_piece_processor = spm.SentencePieceProcessor()
    sentence_piece_processor.load(model_path)

    with open(input_fasta) as input_handler:
        for record in SeqIO.parse(input_handler, format="fasta"):
            tokenized_sequence = sentence_piece_processor.EncodeAsIds(str(record.seq))
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
