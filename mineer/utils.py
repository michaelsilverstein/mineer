"""
minEER utilites
"""
from typing import List, Tuple
from .mineer import minEER
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
import os, functools, random

default_nreads = 5000
default_mal = 100
default_mae = 1e-2

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

        # minEER methods
        self.mineer = False
        self.pass_qc_mineer = None
        self.trimpos_mineer = None

        # Truncation methods
        self.trimmed = None
        self.pass_qc = None
        self.bothpassing = None # Updated from ReadPair

    def runMineer(self):
        """Run minEER on Read"""
        # Run minEER to get truncation positions
        trimpos_mineer = minEER(self.untrimmed.ee, self.mal, self.mae)
        self.trimpos_mineer = trimpos_mineer
        self.mineer = True
        self.pass_qc_mineer = not trimpos_mineer[0] is None

    def truncate(self, trimpos):
        """Truncate sequence given """
        trimstart, trimend = trimpos
        trimlen = trimend - trimstart
        # Truncate
        trimmed_seqrecord = self.untrimmed.record[trimstart: trimend]
        self.trimmed = Record(trimmed_seqrecord)
        # Passes QC?
        self.pass_qc = (self.trimmed.ee.mean() <= self.mae) & (self.trimmed.length == trimlen)


class ReadPair:
    """A forward/rev read pair"""
    def __init__(self, fwd_read: Read, rev_read:Read = None):
        self.fwd_read = fwd_read
        self.rev_read = rev_read
        self.paired = bool(rev_read)

    def truncBoth(self, fwd_pos: tuple, rev_pos: tuple=None):
        """Truncate both read pairs"""
        self.fwd_read.truncate(fwd_pos)
        if self.paired:
            self.rev_read.truncate(rev_pos)
    
    def checkQC(self):
        """Check that both reads are passing qc"""
        # Check that read has been truncated, don't evaluate qc if not
        if not self.fwd_read.trimmed:
            return
        qcs = [self.fwd_read.pass_qc]
        if self.paired:
            if not self.rev_read.trimmed:
                return
            qcs.append(self.rev_read.pass_qc)
        self.bothpassing = all(qcs)
    
    @property
    def bothpassing(self):
        """Check that both reads are passing qc"""
        # Check that read has been truncated, don't evaluate qc if not
        if not self.fwd_read.trimmed:
            return
        qcs = [self.fwd_read.pass_qc]
        if self.paired:
            if not self.rev_read.trimmed:
                return
            qcs.append(self.rev_read.pass_qc)
        
        # Update method for each Read and for ReadPair
        bp = all(qcs)
        self.fwd_read.bothpassing = bp
        if self.paired:
            self.rev_read.bothpassing = bp
        return bp

class Sample:
    """Readpairs for all reads in both files of a sample (for paired reads). Assumes reads are in order for each pair"""
    def __init__(self, name: str, fwd_file: File, rev_file: File=None):
        self.name = name
        self.fwd_file = fwd_file
        if rev_file:
            self.rev_file = rev_file
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
    Collects samples for project

    Inputs:
    | filepaths: Array of file path locations. Files without a file suffix will be ignored
    | {fwd, rev}_format: Suffix of forward and reverse files
        ex. For forward read A_1.fastq and reverse read A_2.fastq
        fwd_format = '_1.fastq'
        rev_format = '_2.fastq'
        * If single ended, do not define rev_format
    | nreads: Number of reads to sample for computing truncation positions (default: 10000)
    | mal: Maximum acceptable length (default: 100)
    | mae: Minimum acceptable error (default: 1e-3)
    | aggmethod: Method to aggregate truncation positions across reads, either 'median' (default) or 'mean'
    | outdir: Output directory

    Pipeline structure
    * Project: Collection of samples
    * Sample: Collection of read pairs
    * ReadPair: Collection of pair of reads
    * Read: Untrimmed and trimmed record
    * Record: A single SeqRecord with extracted data
    """
    def __init__(self, filepaths: list, fwd_format: str, rev_format: str=None, nreads: int=default_nreads, mal: int=default_mal, mae: float=default_mae, aggmethod: str='median', outdir: str=None):
        #TODO: MAKE SURE ABSPATH WORKS
        self.filepaths = [f for f in filepaths if any(map(lambda x: f.endswith(x), [fwd_format, rev_format]))]
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

        # All reads
        self.fwd_reads = []
        self.rev_reads = []
        self.samples = None
        # Subset for minEER
        self.fwd_sub = None
        self.fwd_n_tries = None
        self.fwd_frac_passing = None
        self.rev_sub = None
        self.rev_n_tries = None
        self.rev_frac_passing = None
        self.all_sub = []
        # Global
        self.fwd_pos = None
        self.fwd_len = None
        self.rev_pos = None
        self.rev_len = None
        self.passing_readpairs = None

    @functools.cached_property
    def files(self) -> List[File]:
        """All files in project"""
        fs = []
        for path in self.filepaths:
            direction = 'f' if path.endswith(self.fwd_format) else 'r'
            suffix = self.fwd_format if direction == 'f' else self.rev_format
            file = File(path, direction, suffix, self.mal, self.mae)
            fs.append(file)
        return fs

    def getReadsandSamples(self):
        """Sort reads into forward and reverse reads and create samples"""
        samples = {}
        for file in self.files:
            # Create sample pairs
            sample = file.sample
            direction = file.direction
            if sample not in samples:
                samples[sample] = {'f': None, 'r': None}
            samples[sample][direction] = file

            # Sort reads by direction
            reads = file.reads
            if direction == 'f':
                attr = 'fwd_reads'
            else:
                attr = 'rev_reads'
            for read in reads:
                # Save read to direction
                getattr(self, attr).append(read)

        # Create sample objects
        Samples = [Sample(sample, files['f'], files['r']) for sample, files in samples.items()]
        self.samples = Samples

    def subsampleReads(self, reads: List[Read]):
        """Subsample `nreads` and keep those that pass minEER"""
        # Shuffle reads
        random.shuffle(reads)

        # Keep track of passing
        read_count = 0
        passing_reads = []

        # Get passing reads until `nreads` is reached
        for read in reads:
            read_count += 1
            read.runMineer()
            if read.pass_qc_mineer:
                passing_reads.append(read)
            
            # Break if `nreads` has been reached
            if len(passing_reads) == self.nreads:
                break
        return read_count, passing_reads
    
    def subsampleAll(self):
        """Subsample all reads up to `nreads` keeping those that pass minEER and recording fraction that pass"""
        self.fwd_n_tries, self.fwd_sub = self.subsampleReads(self.fwd_reads)
        self.fwd_frac_passing = len(self.fwd_sub) / self.fwd_n_tries
        self.all_sub.extend(self.fwd_sub)
        if self.paired:
            self.rev_n_tries, self.rev_sub = self.subsampleReads(self.rev_reads)
            self.rev_frac_passing = len(self.fwd_sub) / self.rev_n_tries
            self.all_sub.extend(self.rev_sub)
    
    def _reportPassingSubset(self):
        report = f'From subset:\n\tForward reads passing: {len(self.fwd_sub)} /{self.fwd_n_tries} ({self.fwd_frac_passing * 100:.2f})%\n'
        if self.paired:
            report += f'\tReverse reads passing: {len(self.rev_sub)} /{self.rev_n_tries} ({self.rev_frac_passing * 100:.2f})%\n'
        print(report)

    def calcPos(self):
        """Calculate global truncation positions from subset"""
        # Get forward trim positions
        fwd_positions = [r.trimpos_mineer for r in self.fwd_sub if r.pass_qc_mineer]
        self.fwd_passing_sub = len(fwd_positions)/self.nreads
        self.fwd_pos = self.aggfunc(fwd_positions, 0).astype(int)
        self.fwd_len = self.fwd_pos[1] - self.fwd_pos[0]
        if self.paired:
            rev_positions = [r.trimpos_mineer for r in self.rev_sub if r.pass_qc_mineer]
            self.rev_passing_sub = len(rev_positions)/self.nreads
            self.rev_pos = self.aggfunc(rev_positions, 0).astype(int)
            self.rev_len = self.rev_pos[1] - self.rev_pos[0]

    def truncAndFilter(self):
        """Truncate all reads to global positions and indicate passing ReadPairs"""
        assert not self.fwd_pos is None, 'Run `calcPos()` for global truncation parameters (or set `fwd_pos` and `rev_pos` manually)'

        passing_rps = []
        for sample in self.samples:
            readpairs = sample.readpairs
            for rp in readpairs:
                rp.truncBoth(self.fwd_pos, self.rev_pos)
                # Save read pairs that both pass QC
                if rp.bothpassing:
                    passing_rps.append(rp)
        
        self.passing_readpairs = passing_rps

    def writeFile(self, file: File):
        """Write one file with truncated sequences. Use File to preserve original order"""
        # Generate output path
        outpath = os.path.join(self.outdir, f'{file.sample}_mineer{file.suffix}')
        # Create output directory if it doesn't exist
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        # Gather truncated reads where pairs pass qc
        truncated = [r.trimmed.record for r in file.reads if r.bothpassing]
        # Save
        SeqIO.write(truncated, outpath, 'fastq')
            
    def writeFiles(self):
        """Write all truncated files"""
        for f in self.files:
            self.writeFile(f)
