# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import sys
import click
import pathlib
import logging

from Bio import SeqIO

from gcgc import __version__
from gcgc.encoded_seq import EncodedSeq
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
@click.argument("organism_id")
@click.argument("data_directory")
def download_organism(organism_id, data_directory):
    path = pathlib.Path(data_directory)
    taxon_dataset = dataset.TaxonDataset([83333], path)
    taxon_dataset.download_files()


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
