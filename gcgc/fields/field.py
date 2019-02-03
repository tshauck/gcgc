# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""The base Field type."""


class Field:
    """The base class for a Field, this holds the name."""

    def __init__(self, name: str) -> None:
        """Init the field object."""
        self.name = name
