"""
minEER pipeline

AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
from .utils import Project, default_nreads, default_mal, default_mae


def truncPipeline(filepaths: list, fwd_format: str, rev_format = None, mal=default_mal, mae=default_mae, aggmethod: str='median', nreads=default_nreads, outdir=None):
    """
    Leave `rev_format` blank for single read mode

    Run minEER truncation pipeline
    1) Ingest files
    2) Run minEER on subset
    3) Determine global truncation positions
    4) Truncate all reads to global positions and filter out read pairs that don't pass QC
    5) Save truncated sequences
    """
    # Create project
    print('Creating project...')
    project = Project(filepaths, fwd_format, rev_format, nreads, mal, mae, aggmethod, outdir)
    # 1) Ingest reads
    print('Ingesting reads...')
    project.getReadsandSamples()
    # 2) Subsample
    print('Running minEER on a subset of %d reads per direction...' % nreads)
    project.subsampleAll()
    # 3) Determine global truncation positions
    print('Calculating global truncation positions...')
    project.calcPos()
    # 4) Truncate all reads to global positions and filter out read pairs that don't pass QC
    print('Truncating files to global positions...')
    project.truncAndFilter()
    # 5) Save truncated sequences
    print('Saving truncated files...')
    project.writeFiles()