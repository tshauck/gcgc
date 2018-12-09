# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import logging
import pathlib
import sys

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
    print(__version__)


@main.command()
@click.argument("organism_id", nargs=-1)
@click.argument("data_directory", nargs=1)
def download_organism(organism_id, data_directory):
    path = pathlib.Path(data_directory)
    taxon_dataset = dataset.TaxonDataset(organism_id, path)
    taxon_dataset.download_files()


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
