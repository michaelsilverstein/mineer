"""
minEER utilites
"""
from .mineer import minEER
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
import pandas as pd
import os, functools

default_nreads = 10000
default_mal = 100
default_mae = 1e-3

defaults = dict

class File:
    """Fastq file"""
    def __init__(self, filepath: str, direction: str, suffix: str=None, mal: int=default_mal, mae: float=default_mae):
        self.filepath = filepath
        assert direction in ['f', 'r'], 'Direction must be "f" or "r"'
        self.direction = direction
        self.suffix = suffix
        self.sample = os.path.basename(filepath).split(suffix)[0]
        self.mal = mal
        self.mae = mae

    @functools.cached_property
    def reads(self):
        """Load all reads"""
        return [Read(r, self, self.mal, self.mae) for r in SeqIO.parse(self.filepath, 'fastq')]


class Record:
    """A single read"""
    def __init__(self, seqrecord):
        self.record = seqrecord
        
    @functools.cached_property
    def length(self):
        """Read length"""
        return len(self.record)

    @functools.cached_property
    def phred(self):
        """Extract phred quality scores"""
        return self.record.letter_annotations['phred_quality']
    
    @functools.cached_property
    def ee(self):
        """Get expected error scores"""
        return phred2ee(self.phred)
    
class Read:
    """Untrimmed and Trimmed data on read"""
    def __init__(self, seqrecord: SeqRecord, file: File=None, mal: int=default_mal, mae: float=default_mae):
        self.untrimmed = Record(seqrecord)
        self.id = seqrecord.id
        self.file = file
        self.mal = mal
        self.mae = mae
        self.mineer = False
        self.pass_qc = None
        self.trimmed = None
    
    @functools.cached_property
    def trimpos(self):
        """Get trim positions using minEER"""
        self.mineer = True
        return minEER(self.untrimmed.ee, self.mal, self.mae)
    
    def trim(self, trimpos: tuple=None):
        """Trim read using minEER"""
        # Get trim positions using minEER if none provided
        if trimpos is None:
            trimpos = self.trimpos
        trimstart, trimend = trimpos

        # Check if QC is passed (if minEER was run)
        if (trimstart is None) & (trimend is None):
            self.pass_qc = False
        else:
            # Trim sequence
            trimstart, trimend = map(int, (trimstart, trimend))
            trimmed_seqrecord = self.untrimmed.record[trimstart: trimend + 1]
            self.trimmed = Record(trimmed_seqrecord)
            # Evaluate pass_qc if mineer was not run
            if self.mineer:
                self.pass_qc = True
            else:
                self.pass_qc = self.trimmed.ee.mean() <= self.mae

class ReadPair:
    """A forward/rev read pair"""
    def __init__(self, fwd_read: Read, rev_read:Read = None):
        self.fwd_read = fwd_read
        self.rev_read = rev_read
    
    @functools.cached_property
    def bothpassing(self):
        """Check that both reads are passing"""
        # If not rev read
        if self.rev_read:
            check = self.fwd_read.pass_qc & self.rev_read.pass_qc
        else:
            check = self.fwd_read.pass_qc
        return check

class Sample:
    """A (paired) sample"""
    def __init__(self, name: str, fwd_file: File, rev_file: File=None):
        self.name = name
        self.fwd_file = fwd_file
        if rev_file:
            self.rev_file = rev_file
            # Collect readpairs for 
            self.readpairs = [ReadPair(f, r) for f, r in zip(self.fwd_file.reads, self.rev_file.reads)]
        else:
            self.readpairs = [ReadPair(f) for f in self.fwd_file.reads]
def phred2ee(phred):
    """
    Convert a Phred score to an expected error probability
    ee = 10 ^ (-P/10)
    """
    phred = np.array(phred)
    ee = np.power(10, -phred / 10)
    return ee

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

    def truncAll(self):
        """
        Truncate all reads using global positions
        """
        all_reads = self.fwd_reads.extend(self.rev_reads)
        for r in all_reads:
            truncpos = tuple(map(int, self.fwd_pos if r.file.direction == 'f' else self.rev_pos))
            if not r.mineer:
                r.trim(truncpos)
            # If mineer has been run, but different position, re-trim
            elif r.trimpos != truncpos:
                r.trim(truncpos)

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
    
    def writeSample(self, sample: Sample):
        """Write one sample, ensuring both pairs pass"""
        # Gather passing readpairs
        pass
            
    def writeFiles(self):
        """Write all truncated files"""
        for f in self.files:
            self.writeFile(f)

# def sample2file(all_filepaths, suffix):
#     """
#     Get mapping of sample -> filepath given a suffix
#     """
#     # Get filepaths of `suffix`
#     filepaths = [f for f in all_filepaths if f.endswith(suffix)]
#     # Get all forward sample names
#     samples = [os.path.basename(f).split(suffix)[0] for f in filepaths]
#     # Combine
#     return dict(zip(samples, filepaths))
    

# def pairSamples(files: l, fwd_format: str, rev_format: str=None) -> pd.DataFrame:
#     """
#     Given a list of filepaths and suffixes, create filepath pairs
#     Leave `rev_format` undefined for single end

#     Output: | [Sample] | forward | reverse |
#     """
#     # Get mapping of samples to forward reads
#     fwd_files = [f for f in files if f.]
#     fwd_samples = {os.path.basename(f.filepath).split(fwd_format)[0]: f for f in fwd_files}
#     fwd_df = pd.Series(fwd_samples, name='forward').rename_axis('Sample')
#     fwd_samples = {}
#     for f in files:
#         if f.filepath
#     # If paired, get reverse
#     if rev_format:
#         rev_samples = samples2paths(filepaths, rev_format)
#         rev_df = pd.Series(rev_samples, name='reverse').rename_axis('Sample')
#         # Combine
#         file_df = pd.concat((fwd_df, rev_df), 1)
#     else:
#         # Otherwise just have forward
#         file_df = fwd_df.reset_index()
#     return file_df