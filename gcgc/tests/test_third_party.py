# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Tests for third party modules."""


import pytest

from gcgc.third_party import hf_datasets
from gcgc.third_party import hf_tokenizer


def test_tokenizer():
    """Test the tokenizer interface."""
    tokenizer = hf_tokenizer.build_hf_tokenizer(1, 1, "extended_protein")

    # Includes start and stop ids.
    expected = [11, 18, 11]
    encoded = tokenizer.encode("MVM")
    assert expected == encoded.ids


@pytest.mark.skip
def test_download_swissprot():
    """Test it's possible to download the swissprot dataset."""
    ref = hf_datasets.UniprotDataset(name="sprot")
    ref.download_and_prepare()
    dataset = ref.as_dataset()
    assert "sprot" in dataset.keys()
