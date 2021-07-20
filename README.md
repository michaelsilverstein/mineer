# minEER
A sensible sequence trimming algorithm

# Install
After cloning run `pip install .`

# Tutorial
1. Download some fastq files

```bash
fastq-dump SRR9660307 SRR9660321 --split-files -O test_files
```

2. Run  pipeline
```bash
mineer -i test_files -f _1.fastq -r _1.fastq -o test_out
```
# Methods

Method:
1) Subsample reads from each direction
2) Apply minEER to each subsampled read to determine optimal truncation positions
3) Compute global truncation positions
4) Adjust all reads to global positions and filter out reads that do not meet quality criteria

Assumptions:
* Reads in paired files are ordered

# Contributing
Run tests with `python -m unittest`