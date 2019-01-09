# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""Fixtures for testing."""

import pathlib

_PATH = pathlib.Path(__file__).parent

P53_HUMAN = _PATH / "p53_human/p53_human.fasta"
ECOLI_PATH = _PATH / "ecoli"
PF01152_PATH = _PATH / "globin_alignment/PF01152_seed.sto"
PF01152_PATH_FULL = _PATH / "globin_alignment/PF01152_full.fasta"
PF12057_PATH = _PATH / "PF12057/PF12057.fasta"
