# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Sets up the GCGC package."""

from pathlib import Path

from setuptools import setup, find_packages

DESCRIPTION = "GCGC is a preprocessing library for biological sequence model development."
THIS_DIR = Path(__file__).resolve().parent


def _get_long_description():
    history = (THIS_DIR / "CHANGELOG.md").read_text()
    return (THIS_DIR / "README.md").read_text() + "\n\n" + history


setup(
    name="gcgc",
    version="0.12.0-dev.4",
    description=DESCRIPTION,
    long_description=_get_long_description(),
    long_description_content_type="text/markdown",
    classifiers=["License :: OSI Approved :: MIT License"],
    author="Trent Hauck",
    author_email="trent@trenthauck.com",
    url="https://github.com/tshauck/gcgc",
    license="MIT",
    packages=find_packages(),
    install_requires=["pydantic~=1.1"],
    python_requires=">=3.6",
    extras_require={"sentencepiece": ["sentencepiece~=0.1", "biopython"]},
    entry_points="""
        [console_scripts]
        gcgc=gcgc:cli
    """,
)
