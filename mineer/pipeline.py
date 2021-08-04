"""
minEER pipeline

Github: https://github.com/michaelsilverstein/minEER
AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
import argparse
import random, os
from .utils import Project, default_nreads, default_mal, default_mae
from .viz import Viz
from argparse import ArgumentParser, RawTextHelpFormatter
import pandas as pd


def truncPipeline(filepaths: list, fwd_format: str, rev_format = None, mal=default_mal, mae=default_mae, aggmethod: str='median', nreads=default_nreads, filter='both', outdir=None, viz_outdir=None, no_shuffle=False, write=True):
    """
    Leave `rev_format` blank for single read mode

    See mineer.utils.Project for documentation on each argument

    Run minEER truncation pipeline
    1) Ingest files
    2) Run minEER on subset
    3) Determine global truncation positions
    4) Truncate all reads to global positions and filter out read pairs that don't pass QC
    5) Save truncated sequences
    6) Produce viz if directory provided
    """
    print('******** STARTING MINEER PIPELINE ********')
    # Create project
    project = Project(filepaths, fwd_format, rev_format, nreads, mal, mae, aggmethod, filter, outdir, no_shuffle)
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
    
    # Only truncate and save if `write`
    if write:
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
    parser.add_argument('-i', help='Path to directory containing all fastq files (raw or gzipped) for a single project', required=True, dest='filepaths')
    parser.add_argument('-f', help="""Forward read filename suffix for paired fastqs or single read filename suffix
        Ex. sample1_1.fastq
            -f _1.fastq
        If extension contains a '-', like sample1-r1.fastq:
            -f="-r1.fastq"
    """, required=True, dest='fwd_format')
    parser.add_argument('-r', help='Reverse read filename suffix (leave blank for single end)', dest='rev_format')
    parser.add_argument('--mal', help='Minimum acceptable length. Default: %d' % default_mal, type=int, default=default_mal)
    parser.add_argument('--mae', help='Maximum acceptable error. Deafult: %f' % default_mae, type=float, default=default_mae)
    parser.add_argument('-m', help='Aggregation method for computing truncation. Default: "median"', choices=['mean', 'median'], default='median', dest='aggmethod')
    parser.add_argument('-n', help='Number of reads to subsample per direction for computing truncation position. Default: %d' % default_nreads, type=int, default=default_nreads, dest='nreads')
    parser.add_argument('--filter', help="""How to filter out truncated reads
        --filter any: Filter out read pairs where any read has EER > mae
        --filter both [Default]: Only filter read pairs where both reads have EER > mae
        --filter no: Do not filter read pairs based on EER
        * In all cases, reads that fall within truncation positions will be filtered""", choices=['any', 'both', 'no'], default='both')
    parser.add_argument('-o', help='Output directory. Default: current working directory', dest='outdir')
    parser.add_argument('-v', help='Provide output directory to generate and visualizations', dest='viz_outdir')
    parser.add_argument('--test', help=argparse.SUPPRESS, action='store_true', default=False, dest='no_shuffle')

    args = parser.parse_args(args)
    
    # Get filepaths
    filepaths = [os.path.abspath(os.path.join(args.filepaths, f)) for f in os.listdir(args.filepaths)]

    # Extract inputs
    inputs = args.__dict__
    inputs['filepaths'] = filepaths
    inputs['write'] = bool(args.outdir)

    # Run pipeline
    project = truncPipeline(**inputs)

    return project

def run():
    # Seperate function to access `project` when testing, but avoid returning object to command line
    mineer_cli()
