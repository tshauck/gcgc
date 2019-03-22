# (c) Copyright 2018 Trent Hauck
# All Rights Reserved
"""A set of Datasets which are comprised of DataFiles."""

import asyncio
import logging
from pathlib import Path
from typing import List, Optional, Sequence

import aiohttp

logger = logging.getLogger(__name__)


class DataFile:
    """An external file download and place at the local path."""

    def __init__(self, url: Optional[str], local_path: Path) -> None:
        """Init the DataFile object."""

        self.url = url
        self.local_path = local_path

    async def download(self, session):
        """Given a session, query the url."""

        async with session.get(self.url) as response:
            with open(self.local_path, "wb") as f_handle:
                while True:
                    chunk = await response.content.read(1024)
                    if not chunk:
                        break
                    f_handle.write(chunk)

            logger.info("Downloaded %s to %s.local_path}.", self.url, self.local_path)
            return await response.release()


DataFiles = List[DataFile]


class UniprotDataset:
    """A base class that loads data from Uniprot."""

    def __init__(self, data_files: DataFiles):
        """Initialize the UniprotDataset object."""

        self.data_files = data_files

    @staticmethod
    def _url(query_str, file_format):
        return f"https://www.uniprot.org/uniprot/?{query_str}&format={file_format}"

    def download_files(self):
        """Download the files."""
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait([self._download()]))
        loop.close()

    async def _download(self):
        async with aiohttp.ClientSession() as session:
            downloaded_files = [f.download(session) for f in self.data_files]
            await asyncio.gather(*downloaded_files)


class TaxonDataset(UniprotDataset):
    """Download a set of taxons IDs from uniprot."""

    def __init__(
        self, taxon_ids: Sequence[int], data_directory: Path = Path("."), file_format: str = "fasta"
    ):
        """Initialize the TaxonDataset object."""

        self.taxon_ids = taxon_ids
        files_to_download: DataFiles = []

        for taxon_id in self.taxon_ids:

            url = self._url(f"query=reviewed:yes+AND+organism:{taxon_id}", file_format=file_format)
            local_path = data_directory / str(f"{taxon_id}.{file_format}")

            files_to_download.append(DataFile(url, local_path))

        super().__init__(files_to_download)
