# minEER
A sensible sequence trimming algorithm

Method:
1) Subsample reads from each direction
2) Apply minEER to each subsampled read to determine optimal truncation positions
3) Compute global truncation positions
4) Adjust all reads to global positions and filter out reads that do not meet quality criteria

Assumptions:
* Reads in paired files are ordered