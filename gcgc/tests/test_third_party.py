# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Tests for third party modules."""

import pathlib
import tempfile
from pathlib import Path

import pytest
from transformers import CONFIG_MAPPING
from transformers import AutoModelWithLMHead
from transformers import DataCollatorForLanguageModeling
from transformers import Trainer
from transformers import TrainingArguments

from gcgc import KmerTokenizer
from gcgc.tests.fixtures import PF01152_PATH_FULL
from gcgc.third_party import hf
from gcgc.third_party import hf_datasets


def test_transformers_model(tmpdir):
    """Test a transformers model can be fit with a GCGC."""
    tmpdir = Path(tmpdir)

    kmer_tokenizer = KmerTokenizer(alphabet="extended_protein")
    transformers_tokenizer = hf.GCGCTransformersTokenizer.from_kmer_tokenizer(kmer_tokenizer)

    dataset = hf.GenomicDataset.from_path(PF01152_PATH_FULL, transformers_tokenizer)
    data_collator = DataCollatorForLanguageModeling(tokenizer=transformers_tokenizer, mlm=True)

    config = CONFIG_MAPPING["bert"](
        vocab_size=len(kmer_tokenizer.vocab),
        num_hidden_layers=2,
        hidden_size=12,
        intermediate_size=12,
        num_attention_heads=2,
        type_vocab_size=1,
    )
    model = AutoModelWithLMHead.from_config(config)

    training_args = TrainingArguments(str(tmpdir), max_steps=10)

    trainer = Trainer(
        model=model,
        args=training_args,
        data_collator=data_collator,
        train_dataset=dataset,
        prediction_loss_only=True,
    )

    trainer.train()


def test_tokenizer():
    """Test the tokenizer interface."""
    kmer_tokenizer = KmerTokenizer(alphabet="extended_protein", conform_length=9)

    transformers_tokenizer = hf.GCGCTransformersTokenizer.from_kmer_tokenizer(kmer_tokenizer)

    with tempfile.TemporaryDirectory() as tempdir:
        temppath = pathlib.Path(tempdir)
        model_path = temppath / "model"
        model_path.mkdir(exist_ok=True)

        transformers_tokenizer.save_pretrained(str(model_path))
        transformers_tokenizer = hf.GCGCTransformersTokenizer.from_pretrained(str(model_path))

    # Includes start and stop ids.
    expected = [1, 15, 22, 15, 2]
    assert expected == transformers_tokenizer.encode("MVM")

    # Trim off the start and stop id.
    assert expected[1:-1] == transformers_tokenizer.convert_tokens_to_ids(list("MVM"))

    assert {
        "input_ids": [1, 15, 22, 15, 2, 0, 0, 0, 0],
        "token_type_ids": [0, 0, 0, 0, 0, 0, 0, 0, 0],
        "attention_mask": [1, 1, 1, 1, 1, 0, 0, 0, 0],
    } == transformers_tokenizer.encode_plus(
        "MVM", max_length=9, pad_to_max_length=True, truncation=True
    )

    # Test decode
    assert transformers_tokenizer.decode([1, 15, 22, 3, 2, 0]) == ">MV#<|"


@pytest.mark.skip
def test_download_swissprot():
    """Test it's possible to download the swissprot dataset."""
    ref = hf_datasets.UniprotDataset(name="sprot")
    ref.download_and_prepare()
    dataset = ref.as_dataset()
    assert "sprot" in dataset.keys()
