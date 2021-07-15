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
        self.length
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