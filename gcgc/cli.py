# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import sys
import click

from gcgc import __version__


@click.group()
def main(args=None):
    """Console script for gcgc."""
    return 0


@main.command("version")
def version():
    print(__version__)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
