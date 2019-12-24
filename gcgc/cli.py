# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""CLI Functions for GCGC."""

import json
import logging
import pathlib
import sys
import shutil
import tempfile

import click

try:
    import sentencepiece as spm
    from Bio import SeqIO

    has_spm = True
except ImportError:
    has_spm = False

from gcgc import __version__

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
@click.option("--vocab_size", type=int, default=200, help="The size of the vocabulary to target.")
@click.option("--model_type", type=str, default="unigram", help="Which model, defaults to unigram.")
@click.option("--max_sequence_length", type=int, default=4192, help="How long of sequence to use.")
def sentencepiece_train(
    input_fasta: str, model_prefix: str, vocab_size: int, model_type: str, max_sequence_length: int
):
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

    if not has_spm:
        raise RuntimeError(f"Do not have the SentencePiece python package, see setup.py.")

    if shutil.which("spm_train") is None:
        raise RuntimeError(f"spm_train is not on the path.")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = pathlib.Path(tmpdir)

        text_file_path = tmppath / "input_textfiles.txt"
        with text_file_path.open("w") as text_lines, input_path.open("r") as input_handler:
            for record in SeqIO.parse(input_handler, "fasta"):
                text_lines.write(f"{str(record.seq)}\n")

        arg_string = " ".join(
            [
                f"--input={str(text_file_path)}",
                f"--model_prefix={model_prefix}",
                f"--vocab_size={vocab_size}",
                f"--model_type={model_type}",
                f"--max_sentence_length={max_sequence_length}",
            ]
        )
        spm.SentencePieceTrainer.Train(arg_string)

        output_vocab = input_path.parent / f"{model_prefix}.vocab"
        output_model = input_path.parent / f"{model_prefix}.model"

        LOG.info("Copying %s to %s", output_vocab.name, str(output_vocab))
        shutil.move(output_vocab.name, str(output_vocab))
        shutil.move(output_model.name, str(output_model))
