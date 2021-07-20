"""
minEER pipeline

Github: https://github.com/michaelsilverstein/minEER
AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
import random, os
from .utils import Project, default_nreads, default_mal, default_mae
from .viz import Viz
from argparse import ArgumentParser, RawTextHelpFormatter
import pandas as pd


def truncPipeline(filepaths: list, fwd_format: str, rev_format = None, mal=default_mal, mae=default_mae, aggmethod: str='median', nreads=default_nreads, outdir=None, viz_outdir=None, random_seed=None):
    """
    Leave `rev_format` blank for single read mode

    Run minEER truncation pipeline
    1) Ingest files
    2) Run minEER on subset
    3) Determine global truncation positions
    4) Truncate all reads to global positions and filter out read pairs that don't pass QC
    5) Save truncated sequences
    6) Produce viz if directory provided
    """
    print('******** STARTING MINEER PIPELINE ********')
    if random_seed:
        random.seed(random_seed)
    # Create project
    project = Project(filepaths, fwd_format, rev_format, nreads, mal, mae, aggmethod, outdir)
    # 1) Ingest reads
    project.getReadsandSamples()
    # Print inputs
    print('\t\t**** INPUTS ****')
    input_report, sample_map_df = project._reportInputs
    print(input_report)
    print('SAMPLE-FILE PAIRS')
    with pd.option_context('display.max_rows', None):
        print(sample_map_df)
    print()

    # 2) Subsample
    print('\t\t**** RUNNING MINEER ALGORITHM ****')
    print('Running minEER on a subset of %d reads per direction...' % nreads)
    project.subsampleAll()
    print('Complete.')
    print()
    # Print truncstats
    print('\t\t**** minEER truncation stats from subset ****')
    passing_stats, mineer_stats = project._reportMineerStats
    print(passing_stats)
    print(mineer_stats.T)
    print()

    # 3) Determine global truncation positions
    project.calcPos()
    print('\t\t**** GLOBAL TRUNCATION POSITIONS ****')
    globpos_report = project._reportGlobalPos
    print(globpos_report)
    print()

    # 4) Truncate all reads to global positions and filter out read pairs that don't pass QC
    print('\t\t**** TRUNCATING ALL FILES TO GLOBAL POSITIONS ****')
    project.truncAndFilter()
    truncstats = project._reportTruncStats
    print(truncstats)
    print()

    # 5) Save truncated sequences
    print('\t\t**** SAVING TRUNCATED FILES ****')
    print(f'Saving files to "{project.outdir}" with format <sample>_mineer<suffix>')
    project.writeFiles()

    # 6) Produce viz if directory provided
    if viz_outdir:
        print('\t\t**** GENERATING FIGURES ****')
        v = Viz(project, viz_outdir)
        v.genFigs()
        print(f'Figures have been saved to "{viz_outdir}"')
    print('******** MINEER RUN COMPLETE ********')

    return project

def mineer_cli(args=None):
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
    parser.add_argument('--mae', help='Maximum acceptable error. Deafult: %f' % default_mae, type=float, default=default_mae)
    parser.add_argument('-m', help='Aggregation method for computing truncation. Default: "median"', choices=['mean', 'median'], default='median')
    parser.add_argument('-n', help='Number of reads to subsample per direction for computing truncation position. Default: %d' % default_nreads, type=int, default=default_nreads)
    parser.add_argument('-o', help='Output directory. Default: current working directory', default=os.getcwd())
    parser.add_argument('-v', help='Provide output directory to generate and visualizations')

    args = parser.parse_args(args)
    
    # Get filepaths
    filepaths = [os.path.abspath(os.path.join(args.i, f)) for f in os.listdir(args.i)]

    # Run pipeline
    project = truncPipeline(filepaths, args.f, args.r, args.mal, args.mae, args.m, args.n, args.o, args.v)