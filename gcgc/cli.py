# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""CLI Functions for GCGC."""

import logging

import click

from gcgc import __version__

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@click.group()
def main():
    """Console script for gcgc."""


@main.command()
def version():
    """Print the version and exit."""
    print(__version__)
