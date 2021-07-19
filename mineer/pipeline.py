"""
minEER pipeline

Github: https://github.com/michaelsilverstein/minEER
AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
import random, os
from .utils import Project, default_nreads, default_mal, default_mae
from argparse import ArgumentParser, RawTextHelpFormatter


def truncPipeline(filepaths: list, fwd_format: str, rev_format = None, mal=default_mal, mae=default_mae, aggmethod: str='median', nreads=default_nreads, outdir=None, random_seed=None):
    """
    Leave `rev_format` blank for single read mode

    Run minEER truncation pipeline
    1) Ingest files
    2) Run minEER on subset
    3) Determine global truncation positions
    4) Truncate all reads to global positions and filter out read pairs that don't pass QC
    5) Save truncated sequences
    """
    if random_seed:
        random.seed(random_seed)
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
    return project

def mineer_cli():
    """Command line interface for running minEER"""
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter, description=__doc__)
    parser.add_argument('-i', help='Path to directory containing all (unzipped!) fastq files for a single project', required=True)
    parser.add_argument('-f', help="""Forward read filename suffix for paired fastqs or single read filename suffix
        Ex. sample1_1.fastq
            -f _1.fastq
        If extension contains a '-', like sample1-r1.fastq:
            -f="-r1.fastq"
    """, required=True)
    parser.add_argument('-r', help='Reverse read filename suffix (leave blank for single end)')
    parser.add_argument('--mal', help='Minimum acceptable length. Default: %d' % default_mal, type=int, default=default_mal)
    parser.add_argument('--mae', help='Minimum acceptable error. Deafult: %f' % default_mae, type=float, default=default_mae)
    parser.add_argument('-m', help='Aggregation method for computing truncation. Default: "median"', choices=['mean', 'median'], default='median')
    parser.add_argument('-n', help='Number of reads to subsample per direction for computing truncation position. Default: %d' % default_nreads, type=int, default=default_nreads)
    parser.add_argument('-o', help='Output directory. Default: current working directory', default=os.getcwd())

    args = parser.parse_args()