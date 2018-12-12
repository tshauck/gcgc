!!! Caution
    GCGC is in active development.

# Homepage

GCGC is a package written in Python for pre-processing for biological sequences. Think of it like a
Natural Language Processing pre-processing toolkit with design choices oriented towards the
differences the sequences found in natural language vs biology.

GCGC has two main points of entry. First, imported as a Python package, and,
second, used as a command line tool for specific types of transformations. In this case, the ideas
is to integrate GCGC into a larger data processing pipeline.

## Installation

Install GCGC via pip:

```bash
$ pip install gcgc
```

If you'd like to use one of the third party tools, install the related "extras".

```bash
$ pip install gcgc[torch]
```

## Helpful Links

- __Getting Started__ For a full example of using GCGC with a classification model, see the [splice site
  example](./examples/splice-site.md).

- __Bugs or Help__ Please [file an issue](https://github.com/tshauck/gcgc/issues) if you're running into issues for
  some reason.

- __Development Roadmap__ The GCGC development board is hosted on [notion](https://www.notion.so/3649815c53324f01ae03abc99707dc68?v=98d8b29c39544dca9cde8ddc0dd8c98b).

- __Source Code__ [GitHub Repo](https://github.com/tshauck/gcgc)

## Documentation Version

The documentation you're reading was build for version: 0.3.2-dev.
