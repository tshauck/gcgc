# (c) Copyright 2018 Trent Hauck
# All Rights Reserved


import unittest
from click.testing import CliRunner
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from gcgc.encoded_seq import EncodedSeq
from gcgc.alphabet.iupac import ExtendedIUPACDNAEncoding
from gcgc import cli
from gcgc.tests.fixtures import P53_HUMAN


class TestGcgc(unittest.TestCase):
    """Tests for `gcgc` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        help_result = runner.invoke(cli.main, ["--help"])
        assert help_result.exit_code == 0
        assert "--help  Show this message and exit." in help_result.output

    @unittest.skip("Test failing as invoke isn't writing the file for some reason.")
    def test_convert_to_tf_records(self):

        runner = CliRunner()
        fmt = P53_HUMAN.suffix.lstrip(".")

        runner.invoke(cli.convert_file_to_tf_records, str(P53_HUMAN), fmt)
