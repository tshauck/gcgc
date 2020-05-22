# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Tests for third party modules."""

from pathlib import Path

from transformers import CONFIG_MAPPING
from transformers import AutoModelWithLMHead
from transformers import DataCollatorForLanguageModeling
from transformers import Trainer
from transformers import TrainingArguments

from gcgc import KmerTokenizer
from gcgc import third_party
from gcgc.tests.fixtures import PF01152_PATH_FULL


def test_transformers_model(tmpdir):
    """Test a transformers model can be fit with a GCGC."""
    tmpdir = Path(tmpdir)

    kmer_tokenizer = KmerTokenizer(alphabet="extended_protein")
    transformers_tokenizer = third_party.GCGCTransformersTokenizer(kmer_tokenizer)

    dataset = third_party.GenomicDataset.from_path(PF01152_PATH_FULL, transformers_tokenizer)
    data_collator = DataCollatorForLanguageModeling(tokenizer=transformers_tokenizer, mlm=True,)

    config = CONFIG_MAPPING["bert"](
        vocab_size=len(kmer_tokenizer.vocab),
        num_hidden_layers=2,
        hidden_size=12,
        intermediate_size=12,
        num_attention_heads=2,
        type_vocab_size=1,
    )
    model = AutoModelWithLMHead.from_config(config)

    training_args = TrainingArguments(str(tmpdir))

    trainer = Trainer(
        model=model,
        args=training_args,
        data_collator=data_collator,
        train_dataset=dataset,
        prediction_loss_only=True,
    )

    trainer.train()
