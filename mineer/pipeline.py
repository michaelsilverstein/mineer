"""
minEER pipeline

AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
from .utils import Project, default_nreads, default_mal, default_mae


def truncPipeline(filepaths: list, fwd_format: str, rev_format = None, mal=default_mal, mae=default_mae, aggmethod: str='median', nreads=default_nreads, outdir=None):
    """
    Run minEER truncation pipeline
    1) Ingest files
    2) Subsample to `nreads` per direction
    3) Run minEER on subset
    4) Determine global truncation positions
    5) Truncate all reads to global positions and filter out read pairs that don't pass QC
    6) Save truncated sequences
    """
    # Create project
    print('Creating project...')
    project = Project(filepaths, fwd_format, rev_format, nreads, mal, mae, outdir='truncated')
    # 1) Ingest reads
    print('Ingesting reads...')
    project.getReadsandSamples()
    # 2) Subsample
    print('Subsetting reads...')
    project.sampleReads()
    read_subset = project.all_sub
    # 3) Run minEER on subset
    print('Running minEER on subset...')
    project.mineerReads(read_subset)
    # 4) Determine global truncation positions
    print('Calculating global truncation positions...')
    project.calcPos()
    # 5) Truncate all reads to global positions and filter out read pairs that don't pass QC
    print('Truncating files to global positions...')
    project.truncAndFilter()
    # # 6) Save truncated sequences
    print('Saving truncated files...')
    project.writeFiles()