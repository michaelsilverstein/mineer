"""
minEER pipeline

AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
import random, functools

from Bio import SeqIO
from .utils import File, Sample, default_nreads, default_mal, default_mae
import pandas as pd
import numpy as np
import os

class Project:
    """
    Collects files for project

    Inputs:
    | filepaths: Array of file path locations. Files without a file suffix will be ignored
    | {fwd, rev}_format: Suffix of forward and reverse files
        ex. For forward read A_1.fastq and reverse read A_2.fastq
        fwd_format = '_1.fastq'
        rev_format = '_2.fastq'
        * If single ended, do not define rev_format
    | nreads: Number of reads to sample for computing truncation positions (default: %d)
    | mal: Maximum acceptable length (default: %d)
    | mae: Minimum acceptable error (default: %f)
    | aggmethod: Method to aggregate truncation positions across reads, either 'median' (default) or 'mean'
    | outdir: Output directory
    """ % (default_nreads, default_mal, default_mae)
    def __init__(self, filepaths: list, fwd_format: str, rev_format: str=None, nreads: int=default_nreads, mal: int=default_mal, mae: float=default_mae, aggmethod: str='median', outdir: str=None):
        self.filepaths = [os.path.abspath(f) for f in filepaths if any(map(lambda x: f.endswith(x), [fwd_format, rev_format]))]
        self.fwd_format = fwd_format
        self.rev_format = rev_format
        self.nreads = nreads
        self.paired = bool(rev_format)
        self.mal = mal
        self.mae = mae
        assert aggmethod in ['mean', 'median'], '"aggmethod" must be either "median" (default) or "mean"'
        self.aggfunc = np.median if aggmethod == 'median' else np.mean
        if not outdir:
            outdir = os.getcwd()
        self.outdir = outdir

        self.fwd_reads = []
        self.rev_reads = []
        self.fwd_sub = None
        self.rev_sub = None
        self.all_sub = []
        self.fwd_passing = None
        self.rev_passing = None
        self.fwd_pos = None
        self.rev_pos = None

    @functools.cached_property
    def files(self) -> list:
        """All files in project"""
        fs = []
        for path in self.filepaths:
            direction = 'f' if path.endswith(self.fwd_format) else 'r'
            suffix = self.fwd_format if direction == 'f' else self.rev_format
            file = File(path, direction, suffix, self.mal, self.mae)
            fs.append(file)
        return fs

    @functools.cached_property
    def samples(self) -> dict:
        """Get sample pairs"""
        # Create sample pairs
        samples = {}
        for f in self.files:
            filepath = f.filepath
            if filepath.endswith(self.fwd_format):
                suffix = self.fwd_format
                direction = 'f'
            elif filepath.endswith(self.rev_format):
                suffix = self.rev_format
                direction = 'r'
            # Get sample name
            sample = os.path.basename(filepath).split(suffix)[0]
            if sample not in samples:
                samples[sample] = {'f': None, 'r': None}
            samples[sample][direction] = f
        
        # Create sample objects
        Samples = [Sample(sample, files['f'], files['r']) for sample, files in samples.items()]
        return Samples

    def getReads(self):
        """Sort reads into forward and reverse reads"""
        for file in self.files:
            reads = file.reads
            if file.direction == 'f':
                attr = 'fwd_reads'
            else:
                attr = 'rev_reads'
            for read in reads:
                getattr(self, attr).append(read)

    def sampleReads(self):
        """Subsample `nreads` from (each) direction"""
        self.fwd_sub = random.sample(self.fwd_reads, self.nreads)
        self.all_sub.extend(self.fwd_sub)
        if self.paired:
            self.rev_sub = random.sample(self.rev_reads, self.nreads)
            self.all_sub.extend(self.rev_sub)
    
    def trimSubset(self):
        """Trim all subsampled reads"""
        for r in self.all_sub:
            r.trim()
    
    def calcPos(self):
        """Calculate global truncation positions"""
        # Get forward trim positions
        fwd_positions = [r.trimpos for r in self.fwd_sub if r.pass_qc]
        self.fwd_passing = len(fwd_positions)/self.nreads
        self.fwd_pos = self.aggfunc(fwd_positions, 0)
        if self.paired:
            rev_positions = [r.trimpos for r in self.rev_sub if r.pass_qc]
            self.rev_passing = len(rev_positions)/self.nreads
            self.rev_pos = self.aggfunc(rev_positions, 0)

    def writeFile(self, file: File):
        """Write one file with truncated sequences"""
        # Generate output path
        outpath = os.path.join(self.outdir, f'{file.sample}{file.suffix}')
        # Create output directory if it doesn't exist
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        # Gather truncated reads
        truncated = []
        truncpos = tuple(map(int, self.fwd_pos if file.direction == 'f' else self.rev_pos))
        for r in file.reads:
            # Truncate if mineer has not been run
            if not r.mineer:
                r.trim(truncpos)
            # If mineer has been run, but different position, re-trim
            elif r.trimpos != truncpos:
                r.trim(truncpos)
            # If read passes QC add to truncated list
            if r.pass_qc:
                truncated.append(r.trimmed.record)
        # Save
        SeqIO.write(truncated, outpath, 'fastq')
            
    def writeFiles(self):
        """Write all truncated files"""
        for f in self.files:
            self.writeFile(f)

def main(filepaths: list, fwd_format: str, rev_format = None, mal=default_mal, mae=default_mae, aggmethod: str='median', nreads=default_nreads, outdir=None):
    """
    Run minEER truncation pipeline
    1) Ingest files
    2) Subsample to `nreads` per direction
    3) Get truncation position for each subsampled read
    4) Determine global truncation positions
    5) Save truncated sequences
    """
    # Create project
    print('Creating project...')
    project = Project(filepaths, fwd_format, rev_format, nreads, mal, mae, aggmethod, outdir)
    # 1) Ingest reads
    print('Ingesting reads...')
    project.getReads()
    # 2) Subsample
    print('Subsetting reads...')
    project.sampleReads()
    # 3) Get truncation position for each subsampled read
    print('Truncating subset...')
    project.trimSubset()
    # 4) Determine global truncation positions
    print('Calculating global truncation positions...')
    project.calcPos()
    # 5) Save truncated sequences
    print('Saving truncated files...')
    project.writeFiles()