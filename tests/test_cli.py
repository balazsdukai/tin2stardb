#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the CLI."""

import pytest

from click.testing import CliRunner

from tin import main


class TestCLI:
    def test_help(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(main.main)
        assert result.exit_code == 0
        assert 'Read/Write a TIN to a database from files.' in result.output
        help_result = runner.invoke(main.main, ['--help'])
        assert help_result.exit_code == 0
        assert 'Usage: main [OPTIONS] CONFIGURATION' in help_result.output