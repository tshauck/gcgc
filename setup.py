# (c) Copyright 2019 Trent Hauck
# All Rights Reserved
"""Sets up the GCGC package."""

from pathlib import Path

from setuptools import setup, find_packages

DESCRIPTION = "GCGC is a preprocessing library for biological sequence model development."
THIS_DIR = Path(__file__).resolve().parent


def get_long_description():
    history = (THIS_DIR / "CHANGELOG.md").read_text()
    return (THIS_DIR / "README.md").read_text() + "\n\n" + history


INSTALL_REQUIRES = [
    "biopython>=1.72",
    "numpy>=1.15.2",
    "click>=7.0",
    "idna_ssl>=1.1",
]

# torch = { version = ">=1.0.0", optional = true }

setup(
    name="gcgc",
    version="0.9.2-dev.1",
    description=DESCRIPTION,
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    classifiers=["License :: OSI Approved :: MIT License"],
    author="Trent Hauck",
    author_email="trent@trenthauck.com",
    url="https://github.com/tshauck/gcgc",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=INSTALL_REQUIRES,
    extras_require={"torch": ["torch>=1.0.0"]},
)
