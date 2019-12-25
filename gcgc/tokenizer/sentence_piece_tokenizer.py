# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Module for Sentence Piece tokenization."""

import tempfile
from typing import List
from pathlib import Path
import shutil

from pydantic import Field
from typing_extensions import Literal

from Bio import SeqIO
from gcgc.tokenizer.base import SequenceTokenizer, SequenceTokenizerSettings

try:
    import sentencepiece as spm

    # pylint: disable=invalid-name
    has_spm = True
except ImportError:
    # pylint: disable=invalid-name
    has_spm = False


class BioSequencePieceSettings(SequenceTokenizerSettings):
    """The settings for the sentence piece model."""

    model_prefix: Path = Field(..., env="GCGC_SP_MODEL_PREFIX")
    vocab_size: int = Field(8000, env="GCGC_SP_VOCAB_SIZE")
    model_type: Literal["unigram", "bpe"] = "unigram"
    max_sequence_length: int = 4192

    @property
    def model_path(self) -> Path:
        """Return the model path based on the prefix."""
        return self.model_prefix.with_suffix(".model")

    @property
    def model_vocab(self) -> Path:
        """Return the model vocab based on the prefix."""
        return self.model_prefix.with_suffix(".vocab")


class BioSequencePiece(SequenceTokenizer):
    """A sentence piece for model on biological sequences."""

    def __init__(self, settings: BioSequencePieceSettings):
        """Init the BioSequencePiece class.

        Args:
            settings: The settings for the tokenizer.

        """
        if not has_spm or not shutil.which("spm_train"):
            raise RuntimeError("Trying to use sentencepiece but the python library is missing!")

        self.settings = settings or BioSequencePieceSettings()

        self.vocab = {}
        self._sp_processor = None

    @property
    def sp_processor(self):
        """Returns the SequencePiece process object."""
        if self._sp_processor is not None:
            return self._sp_processor
        else:
            self._sp_processor = spm.SentencePieceProcessor()
            self._sp_processor.load(str(self.settings.model_path))
            return self._sp_processor

    def fit_on_fasta(self, fasta_file: Path):
        """Run the the SP algo on the fasta_file."""

        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)

            text_file_path = tmppath / "input_textfiles.txt"
            with text_file_path.open("w") as text_lines, fasta_file.open("r") as input_handler:
                for record in SeqIO.parse(input_handler, "fasta"):
                    text_lines.write(f"{str(record.seq)}\n")

            self.fit_on_text(text_file_path)

    def fit_on_text(self, text_file: Path):
        """Run the the SP algo on the text_file."""
        args = [
            f"--input={str(text_file)}",
            f"--model_prefix={self.settings.model_prefix}",
            f"--vocab_size={self.settings.vocab_size}",
            f"--model_type={self.settings.model_type}",
            f"--max_sentence_length={self.settings.max_sequence_length}",
        ]

        vocab_offset = 0
        if self.settings.unk_token:
            self.vocab[vocab_offset] = self.settings.unk_token
            args.extend([f"--unk_piece={self.settings.unk_token}", f"--unk_id={vocab_offset}"])
            vocab_offset += 1
        else:
            args.extend(["--unk_id=-1"])

        if self.settings.bos_token:
            self.vocab[vocab_offset] = self.settings.bos_token
            args.extend([f"--bos_piece={self.settings.bos_token}", f"--bos_id={vocab_offset}"])
            vocab_offset += 1
        else:
            args.extend(["--bos_id=-1"])

        if self.settings.eos_token:
            self.vocab[vocab_offset] = self.settings.eos_token
            args.extend([f"--eos_piece={self.settings.eos_token}", f"--eos_id={vocab_offset}"])
            vocab_offset += 1
        else:
            args.extend(["--eos_id=-1"])

        if self.settings.pad_token:
            self.vocab[vocab_offset] = self.settings.pad_token
            args.extend([f"--pad_piece={self.settings.pad_token}", f"--pad_id={vocab_offset}"])
            vocab_offset += 1
        else:
            args.extend(["--pad_id=-1"])

        spm.SentencePieceTrainer.Train(" ".join(args))

    def encode(self, seq: str) -> List[int]:
        """Encode the underlying sequence into a list of tokens."""
        return self.sp_processor.EncodeAsIds(seq)

    def encode_as_tokens(self, seq: str) -> List[str]:
        """Tokenize the sequence into a list of token tokens.

        Args:
            seq: The sequence to encode.

        Returns:
            The list of strs that are the tokens.

        """
        return self.sp_processor.EncodeAsPieces(seq)

    def load_vocab(self):
        """Load the vocabulary from the file."""
        for line, token in enumerate(self.settings.model_vocab.open()):
            self.vocab[line] = token.strip("\n").split("\t")[0]
