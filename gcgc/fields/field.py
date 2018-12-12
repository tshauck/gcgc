# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""The base Field type."""

from dataclasses import dataclass


@dataclass
class Field:
    """The base class for a Field, this holds the name."""

    name: str
