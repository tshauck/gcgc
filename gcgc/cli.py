# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""CLI Functions for GCGC."""

import logging
import pathlib

import click

from gcgc import __version__
from gcgc.datasets import dataset

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@click.group()
def main(args=None):
    """Console script for gcgc."""
    return 0


@main.command()
def version():
    """Print the version and exit."""
    print(__version__)


@main.command()
@click.argument("organism_id", nargs=-1)
@click.argument("data_directory", nargs=1)
def download_organism(organism_id, data_directory):
    """Download the organisms into the data directory."""

    path = pathlib.Path(data_directory)
    if not path.exists():
        raise Exception(f"The directory {path} does not exist.")

    if not path.is_dir():
        raise Exception(f"Found {path}, but it doesn't appear to be a directory.")

    taxon_dataset = dataset.TaxonDataset(organism_id, path)
    taxon_dataset.download_files()
