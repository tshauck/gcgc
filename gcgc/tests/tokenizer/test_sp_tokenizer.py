# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the sentence models."""

import numpy as np
from gcgc.tokenizer import sentence_piece_tokenizer


def test_train_sentence_piece(tmp_path):
    """Test fitting the sentence piece model."""
    sequences = ["".join(x) for x in np.random.choice(list("ATCG"), (20, 20))]

    text_file = tmp_path / "sp.txt"
    text_file.write_text("\n".join(sequences))

    model_prefix = tmp_path / "model"

    settings = sentence_piece_tokenizer.BioSequencePieceSettings(
        model_prefix=model_prefix, vocab_size=20
    )
    sp_tokenizer = sentence_piece_tokenizer.BioSequencePiece(settings)
    sp_tokenizer.fit_on_text(text_file)

    sp_tokenizer.load_vocab()
    tokenized = sp_tokenizer.encode_as_tokens("ATCGATCGATCG")
    assert all([isinstance(x, str) for x in tokenized])

    ids = sp_tokenizer.encode_as_ids("ATCGATCGATCG")
    assert all([isinstance(x, int) for x in ids])
