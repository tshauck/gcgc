# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import sys
import click
import pathlib
import logging

from Bio import SeqIO
import tensorflow as tf

from gcgc import __version__
from gcgc.encoded_seq import EncodedSeq
from gcgc.third_party.tensorflow_utils import record

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@click.group()
def main(args=None):
    """Console script for gcgc."""
    return 0


@main.command()
def version():
    print(__version__)


def to_path(ctx, param, value) -> pathlib.Path:
    return pathlib.Path(value)


@main.command()
@click.argument("filename", callback=to_path)
@click.argument("format")
def convert_file_to_tf_records(filename, format):

    output_file = filename.with_suffix(".tf_records")
    logger.info(f"Reading from {filename} in format {format} and writing to {output_file}.")

    writer = tf.python_io.TFRecordWriter(str(output_file))
    try:
        with open(filename, "rU") as handle:
            for seq_record in SeqIO.parse(handle, format):
                encoded_seq = EncodedSeq.from_seq(seq_record.seq)
                example = record.to_tensorflow_record(encoded_seq)
                writer.write(example.SerializeToString())
    finally:
        writer.close()


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
