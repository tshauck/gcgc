# (c) Copyright 2018 Trent Hauck
# All Rights Reserved


import unittest

from click.testing import CliRunner

from gcgc import cli


class TestGcgc(unittest.TestCase):
    """Tests for `gcgc` package."""

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        help_result = runner.invoke(cli.main, ["--help"])
        assert help_result.exit_code == 0
        assert "--help  Show this message and exit." in help_result.output
