# (c) Copyright 2020 Trent Hauck
# All Rights Reserved
"""Sequence dataset implementations for the hugging face dataset library.

The hugging face datasets library has objects to help build datasets and work with them. This code
focuses on helping to build datasets from bioinformatics formats such that the rest of the built in
dataset library can be used seamless on bio data.

At present this has uniprot related datasets.

    >>> from gcgc.third_party import hf_datasets
    >>> ref = hf_datasets.UniprotDataset(name="sprot")
    >>> # ref.download_and_prepare()
    >>> # ds = ref.as_dataset()

See the `UniprotDatasetNames` enum for the available names.

    >>> from gcgc.third_party import hf_datasets
    >>> hf_datasets.UniprotDatasetNames.uniref100
    <UniprotDatasetNames.uniref100: 'uniref100'>
"""

# This is redundant, but can be removed when 3.6 support is dropped, and get_args is available.
from enum import Enum
from typing import Dict

try:
    import datasets

    from Bio import SeqIO

except ImportError as exp:
    # pylint: disable=invalid-name
    needed = "datasets"
    raise ImportError(f"Missing one or more libraries: {needed}. Please install: {exp}") from exp


# pylint: disable=line-too-long
# pylint: disable=too-few-public-methods


class FastaBasedBuilder(datasets.GeneratorBasedBuilder):
    """Builds a dataset backed by FASTA files.

    This Builder implements the ability to iterate through a split, but it is incumbant on the
    subclass to implement the _split_generators method.

    When subclassing this, use `self.features` in the dataset info.
    """

    def _split_generators(self, dl_manager: datasets.DownloadManager):
        raise NotImplementedError("Must be implemented by subclass.")

    @property
    def features(self) -> datasets.Features:
        """Return the set of FASTA realted features."""
        return datasets.Features(
            {
                "sequence": datasets.Value("string"),
                "description": datasets.Value("string"),
                "id": datasets.Value("string"),
            }
        )

    def _generate_examples(self, **kwargs):
        """Generate examples from a split.

        This method assumes that kwargs contains file_path. Assuming that's the case, then each
        file_path is looped over, and the contents are yielded as tuples of id a dictionary that
        corresponds to `self.features`.

        Yields:
            Tuples of the id string, and a dictionary with keys "sequence", "description", "id".

        Raises:
            ValueError: In the event the file_paths kwargs isn't passed.

        """
        file_paths = kwargs.get("file_paths")
        if not file_paths:
            raise ValueError("Must pass file_paths.")

        for file_path in file_paths:
            for record in SeqIO.parse(file_path, "fasta"):
                yield record.id, {
                    "sequence": str(record.seq),
                    "description": str(record.description),
                    "id": str(record.id),
                }


_UNIPROT_CITATION = """
@article{10.1093/nar/gky1049,
    author = {The UniProt Consortium},
    title = "{UniProt: a worldwide hub of protein knowledge}",
    journal = {Nucleic Acids Research},
    volume = {47},
    number = {D1},
    pages = {D506-D515},
    year = {2018},
    month = {11},
    issn = {0305-1048},
    doi = {10.1093/nar/gky1049},
    url = {https://doi.org/10.1093/nar/gky1049},
    eprint = {https://academic.oup.com/nar/article-pdf/47/D1/D506/27437297/gky1049.pdf},
}
"""

_UNIPROT_DESCRIPTION = """
The UniProt Knowledgebase is a collection of sequences and annotations for over 120 million proteins across all branches of life. Detailed annotations extracted from the literature by expert curators have been collected for over half a million of these proteins. These annotations are supplemented by annotations provided by rule based automated systems, and those imported from other resources. In this article we describe significant updates that we have made over the last 2 years to the resource. We have greatly expanded the number of Reference Proteomes that we provide and in particular we have focussed on improving the number of viral Reference Proteomes. The UniProt website has been augmented with new data visualizations for the subcellular localization of proteins as well as their structure and interactions. UniProt resources are available under a CC-BY (4.0) license via the web at https://www.uniprot.org/.
    """


class UniprotDatasetNames(str, Enum):
    """Enum for the uniprot datasets."""

    uniref50 = "uniref50"
    uniref90 = "uniref90"
    uniref100 = "uniref100"

    uniparc = "uniparc"
    trembl = "trembl"
    sprot = "sprot"

    def __str__(self) -> str:
        """Use the enum as a string."""
        return f"{self.value}"


class UniprotDatasetConfig(datasets.BuilderConfig):
    """Builder config for uniref datasets."""

    DATASET_URLS: Dict[UniprotDatasetNames, str] = {
        UniprotDatasetNames.uniref50: "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz",
        UniprotDatasetNames.uniref90: "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz",
        UniprotDatasetNames.uniref100: "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz",
        UniprotDatasetNames.sprot: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
        UniprotDatasetNames.trembl: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
        UniprotDatasetNames.uniparc: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz",
    }

    version = datasets.Version("1.0.0")

    def __init__(self, name: UniprotDatasetNames):
        """Init the uniprot dataset with the name of the uniprot dataset."""
        self.name: UniprotDatasetNames = name

    @property
    def url(self) -> str:
        """Return the URL for dataset."""
        return self.DATASET_URLS[self.name]


class UniprotDataset(FastaBasedBuilder):
    """Rpresents a uniprot dataset using the underlying FastaBasedBuilder."""

    VERSION = datasets.Version("1.0.0")

    BUILDER_CONFIG_CLASS = UniprotDatasetConfig
    BUILDER_CONFIGS = [UniprotDatasetConfig(name=name) for name in UniprotDatasetNames]

    def _info(self):
        return datasets.DatasetInfo(
            description=_UNIPROT_DESCRIPTION, features=self.features, citation=_UNIPROT_CITATION
        )

    def _split_generators(self, dl_manager: datasets.DownloadManager):
        """Return the split generators for the uniref datasets."""
        downloaded_file = dl_manager.download_and_extract(self.config.url)
        return [
            datasets.SplitGenerator(
                name=str(self.config.name),
                # These kwargs will be passed to _generate_examples
                gen_kwargs={"file_paths": [downloaded_file]},
            ),
        ]
