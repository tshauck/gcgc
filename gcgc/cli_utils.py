# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Functions that the CLI commands delegate to."""

import pathlib
from typing import Any
from typing import Dict
from typing import Iterator

from Bio import SeqIO

from gcgc.tokenizer.sentence_piece_tokenizer import BioSequencePiece
from gcgc.tokenizer.sentence_piece_tokenizer import BioSequencePieceSettings


def sentencepiece_trainer(input_fasta: str, model_prefix: str):
    """Train the SP model."""
    input_path = pathlib.Path(input_fasta)
    if not input_path.exists():
        raise ValueError(f"{input_path} does not exist.")

    bsp = BioSequencePiece(settings=BioSequencePieceSettings(model_prefix=model_prefix))
    bsp.fit_on_fasta(input_path)


def sentencepiece_tokenizer(input_fasta: str, model_prefix: str) -> Iterator[Dict[str, Any]]:
    """Process the records and yield the results to the CLI.

    Yields:
        Tokenized records.

    """
    tokenizer = BioSequencePiece(settings=BioSequencePieceSettings(model_prefix=model_prefix))

    with open(input_fasta) as input_handler:
        for record in SeqIO.parse(input_handler, format="fasta"):
            tokenized_sequence = tokenizer.encode(str(record.seq))
            output_record = {
                "token_ids": tokenized_sequence,
                "sequence_id": record.id,
                "sequence_description": record.description,
                "n_tokens": len(tokenized_sequence),
            }

            yield output_record
