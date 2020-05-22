# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""The Vocab object holds the vocabulary for the tokenizer."""

import typing

import pydantic

# pylint: disable=unsupported-membership-test
# pylint: disable=unsupported-assignment-operation
# pylint: disable=unsubscriptable-object
# pylint: disable=no-member


class Vocab(pydantic.BaseModel):
    """The Vocab object that holds the mapping between integers and characters."""

    stoi: typing.Dict[str, int] = {}

    @property
    def itos(self):
        """Return the integer to string mapping of the underlying stoi."""
        return {integer: string for string, integer in self.stoi.items()}

    def __getitem__(self, item: str) -> int:
        """Get the item from the vocab."""
        return self.stoi[item]

    def __contains__(self, item: str) -> bool:
        """Return True if the item is in the vocab."""
        return item in self.stoi

    def __len__(self) -> int:
        """Return the size of the vocabulary."""
        return len(self.stoi)

    def __setitem__(self, item: str, value: int) -> None:
        """Add an item to the vocabulary with a specific value."""
        self.stoi[item] = value

    def add_item(self, item: str) -> None:
        """Add an item to the vocabulary."""
        try:
            current_max = max(self.stoi.values())
            self.stoi[item] = current_max + 1
        except ValueError:
            self.stoi[item] = 0

    def add_items(self, items: typing.Iterable[str]) -> None:
        """Add multiple items to the vocabulary in one fell swoop."""
        for item in items:
            self.add_item(item)
