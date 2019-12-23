# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""A vocabulary for tokenization data."""

from typing import Dict

from pydantic import BaseModel


class KmerVocab(BaseModel):
    """A vocabulary object."""

    token_to_int: Dict[str, int]
    int_to_token: Dict[int, str]

    def __len__(self) -> int:
        """Return the length of the vocab."""
        return len(self.token_to_int)
