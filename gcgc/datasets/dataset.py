# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

import asyncio
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

import aiohttp

logger = logging.getLogger(__name__)


@dataclass
class DataFile:
    """
    A an external file download and place at the local path.
    """

    url: Optional[str]
    local_path: Path

    async def download(self, session):
        """
        Given a session, query the url.

        Args:
            session: The asyncio session.

        Returns:
            The released response.
        """

        async with session.get(self.url) as response:
            with open(self.local_path, "wb") as f_handle:
                while True:
                    chunk = await response.content.read(1024)
                    if not chunk:
                        break
                    f_handle.write(chunk)

            logger.info(f"Downloaded {self.url} to {self.local_path}.")
            return await response.release()


DataFiles = List[DataFile]


class UniprotDataset(object):
    """
    A base class that loads data from Uniprot.
    """

    def __init__(self, data_files: DataFiles, file_format: str = "fasta"):
        """
        Initialize the UniprotDataset object.

        Args:
            data_files: The list for FileToDownload object.
            file_format: The file format to use.
        """

        self.data_files = data_files
        self.file_format = file_format

    @staticmethod
    def _url(query_str, file_format):
        return f"https://www.uniprot.org/uniprot/?{query_str}&format={file_format}"

    def download_files(self):
        asyncio.run(self._download())

    async def _download(self):
        async with aiohttp.ClientSession() as session:
            downloaded_files = [f.download(session) for f in self.data_files]
            await asyncio.gather(*downloaded_files)


class TaxonDataset(UniprotDataset):
    """
    Download a set of taxons IDs from uniprot.
    """

    def __init__(
        self, taxon_ids: Sequence[int], data_directory: Path = Path("."), file_format: str = "fasta"
    ):
        """
        Initialize the TaxonDataset object.

        Args:
            taxon_ids: A list of integers representing taxon ids.
            data_directory: The data directory path for the taxon.
            file_format: The file format to use for the downloaded data.
        """

        self.taxon_ids = taxon_ids
        files_to_download: DataFiles = []

        for taxon_id in self.taxon_ids:

            url = self._url(f"query=reviewed:yes+AND+organism:{taxon_id}", file_format=file_format)
            local_path = data_directory / str(f"{taxon_id}.{file_format}")

            files_to_download.append(DataFile(url, local_path))

        super().__init__(files_to_download, file_format)
