# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the sentence models."""

import json
import random
from pathlib import Path

import pytest

from gcgc import cli_utils
from gcgc.alphabets import resolve_alphabet
from gcgc.tests import fixtures
from gcgc.tokenizer import sentence_piece_tokenizer


@pytest.mark.parametrize("alphabet", ["protein", "unambiguous_dna"])
def test_train_sentence_piece(tmp_path, alphabet):
    """Test fitting the sentence piece model."""
    alphabet_list = resolve_alphabet(alphabet, require=True)

    sequences = []
    for _ in range(500):
        new_seq = []
        for _ in range(500):
            new_seq.append(random.choice(list(alphabet_list)))

        sequences.append("".join(new_seq))

    model_prefix = tmp_path / "model"

    settings = sentence_piece_tokenizer.BioSequencePieceSettings(
        model_prefix=model_prefix, vocab_size=50, bos_token=">", eos_token="<", pad_token="|"
    )
    sp_tokenizer = sentence_piece_tokenizer.BioSequencePiece(settings)
    sp_tokenizer.fit_on_list(sequences)

    test_case = "".join(random.choice(alphabet_list) for _ in range(20))

    tokenized = sp_tokenizer.encode_as_tokens(test_case)
    assert all([isinstance(x, str) for x in tokenized])

    ids = sp_tokenizer.encode(test_case)
    assert all([isinstance(x, int) for x in ids])

    assert [sp_tokenizer.vocab[i] for i in tokenized] == ids


def test_cli_functions_e2e(tmp_path):
    """Test we can run on the SP models."""
    prefix = f"{tmp_path}/model"
    cli_utils.sentencepiece_trainer(fixtures.PF01152_PATH_FULL, prefix)
    tokenized_seqs = list(cli_utils.sentencepiece_tokenizer(fixtures.PF01152_PATH_FULL, prefix))
    print(tokenized_seqs)


def test_settings_serialization(tmp_path):
    """Test writing the settings to a file, then getting its results."""
    settings_path = Path(tmp_path) / "settings.json"
    settings = sentence_piece_tokenizer.BioSequencePieceSettings(model_prefix="test")

    settings_path.write_text(json.dumps(settings.dict()))
    settings = json.loads(settings_path.read_text())
