# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Test the alphabet."""

from contextlib import contextmanager

import pytest

from gcgc.tokenizer.base import SequenceTokenizerSettings


@contextmanager
def does_not_raise():
    yield


@pytest.mark.parametrize(
    "settings_dict,expected_token_ids,raises",
    [
        (
            dict(pad_token="|", pad_token_id=1, mask_token="#", mask_token_id=0),
            {"mask_token_id": 0, "pad_token_id": 1},
            does_not_raise(),
        ),
        (
            dict(pad_token="|", mask_token="#", mask_token_id=0),
            {"mask_token_id": 0},
            pytest.raises(ValueError),
        ),
    ],
)
def test_manually_setting_token_ids(settings_dict, expected_token_ids, raises):
    with raises:
        settings = SequenceTokenizerSettings.parse_obj(settings_dict)
        assert settings.mask_token_id == expected_token_ids["mask_token_id"]
        assert settings.pad_token_id == expected_token_ids["pad_token_id"]
