# SAM (BAM) Tools Support

At present, GCGC does not support working with BAM/SAM files directly. Instead one can pre-process
the files by using `samtools fasta` or `samtools fastq`, for example. In particular, see the parsing
options to ensure data from the alignments transfers to the file, then use one of the built-in
GCGC Field parsers.

For example, by running `samtools fasta -T RG ./alignments.bam`, the RG tag from each alignment is
copied into the description of the fasta sequence. Then, by using the DescriptionField in GCGC, the
RG tag can converted into a label. The preprocessing function might look like:

```python
def parse_description(d: str): str:
    for di in d.split("\t"):
        if di.startswith("RG"):
            return di
    else:
        raise RuntimeError("No RG tag found.")
```
