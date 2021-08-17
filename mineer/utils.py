"""
minEER utilites
"""
from typing import List
from .mineer import minEER
import numpy as np
import pandas as pd
import os, functools, random, dnaio
from joblib import Parallel, delayed

default_nreads = 15000
default_mal = 150
default_mae = 1e-2

class File:
    """Fastq file"""
    def __init__(self, filepath: str, direction: str, mal: int, mae: float, suffix: str=None):
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
        rs = []
        with dnaio.open(self.filepath) as fh:
            for r in fh:
                rs.append(Read(r, self.mal, self.mae, self))
        return rs

class Record:
    """A single read"""
    def __init__(self, seqrecord: dnaio._core.Sequence):
        self.record = seqrecord
        
    @functools.cached_property
    def length(self):
        """Read length"""
        return len(self.record)

    @functools.cached_property
    def phred(self):
        """Extract phred quality scores"""
        return [x-33 for x in self.record.qualities.encode()]
    
    @functools.cached_property
    def ee(self):
        """Get expected error scores"""
        return np.array([10 ** (-(x-33)/10)  for x in self.record.qualities.encode()])
    
class Read:
    """Untrimmed and Trimmed data on read"""
    def __init__(self, seqrecord: dnaio._core.Sequence, mal: int, mae: float, file: File=None):
        self.untrimmed = Record(seqrecord)
        self.id = seqrecord.name
        self.file = file
        self.mal = mal
        self.mae = mae

        # minEER methods
        self._trimpos_mineer = None

        # Truncation methods
        self.trimmed = None
        self.pass_e = None
        self.pass_l = None
        self.pass_qc = None
        self.bothpassing = None # Updated from ReadPair
    
    @property
    def trimpos_mineer(self):
        return self._trimpos_mineer
    
    @trimpos_mineer.setter
    def trimpos_mineer(self, trimpos):
        self._trimpos_mineer = trimpos

    @property
    def pass_qc_mineer(self):
        if self.trimpos_mineer is not None:
            return not self.trimpos_mineer[0] is None
        else:
            return None
    
    @property
    def mineer(self):
        """Indicate if mineer has been run"""
        return self.pass_qc_mineer is not None

    def runMineer(self):
        """Run minEER on Read"""
        # Run minEER to get truncation positions
        trimpos_mineer = minEER(self.untrimmed.ee, self.mal, self.mae)
        self._trimpos_mineer = trimpos_mineer

    def truncate(self, trimpos):
        """Truncate sequence given """
        trimstart, trimend = trimpos
        trimlen = trimend - trimstart
        # Truncate
        trimmed_seqrecord = self.untrimmed.record[trimstart: trimend]
        self.trimmed = Record(trimmed_seqrecord)
        # Passes QC?
        self.pass_e = self.trimmed.ee.mean() <= self.mae
        self.pass_l = self.trimmed.length == trimlen
        self.pass_qc = self.pass_e & self.pass_l


class ReadPair:
    """A forward/rev read pair"""
    def __init__(self, fwd_read: Read, rev_read:Read = None, filter: str='both'):
        self.fwd_read = fwd_read
        self.rev_read = rev_read
        self.paired = bool(rev_read)
        self.filter = filter

    def truncBoth(self, fwd_pos: tuple, rev_pos: tuple=None):
        """Truncate both read pairs"""
        self.fwd_read.truncate(fwd_pos)
        if self.paired:
            self.rev_read.truncate(rev_pos)
    
    @property
    def bothpassing(self):
        """Check that both reads are passing qc"""
        # Check that read has been truncated, don't evaluate qc if not
        if not self.fwd_read.trimmed:
            return

        reads = [self.fwd_read]
        if self.paired:
            reads.append(self.rev_read)
        
        # Determine passing using `filter`
        # Filter out readpairs where any read fails. Keep == all passing
        if self.filter == 'any': 
            bp = all([r.pass_qc for r in reads])
        # Filter only when both pass. Keep == any passing 
        if self.filter == 'both':
            bp = any([r.pass_qc for r in reads])
        # Don't filter any based on quality. Keep == all passing length
        if self.filter == 'no':
            bp = all([r.pass_l for r in reads])

        # Update method for each Read and for ReadPair
        self.fwd_read.bothpassing = bp
        if self.paired:
            self.rev_read.bothpassing = bp
        return bp

class Sample:
    """Readpairs for all reads in both files of a sample (for paired reads). Assumes reads are in order for each pair"""
    def __init__(self, name: str, fwd_file: File, rev_file: File=None, filter: str='both'):
        self.name = name
        self.fwd_file = fwd_file
        self.readpairs = None
        if rev_file:
            self.rev_file = rev_file
            self.readpairs = [ReadPair(f, r, filter) for f, r in zip(self.fwd_file.reads, self.rev_file.reads)]
        else:
            self.readpairs = [ReadPair(f, filter=filter) for f in self.fwd_file.reads]
    

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
    | nreads: Minimum number of reads to subsample for computing global truncation positions

    | mal: Maximum acceptable length (default: 100)
    | mae: Maximum acceptable error (default: 1e-2)
    | aggmethod: Method to aggregate truncation positions across reads, either 'median' (default) or 'mean'
    | filter: 
        filter='any':  Filter out all readpairs where any read where EER > mae
        filter='both': Only filter out readpairs where both reads EER > mae [Default]
        filter='no':   Do not filter out any reads based on EER
        * All reads with length < truncation length are always filtered (len(all reads) == truncation length)
    | outdir: Output directory
    | no_shuffle: Don't shuffle reads before sampling (just for testing)
    | nproc: Number of processors to use for parallel operations [Default: -1 (all processors)]

    Pipeline structure
    * Project: Collection of samples
    * Sample: Collection of read pairs
    * ReadPair: Collection of pair of reads
    * Read: Untrimmed and trimmed record
    * Record: A single SeqRecord with extracted data
    """
    def __init__(self, filepaths: List[str], fwd_format: str, rev_format: str=None, nreads: int=default_nreads, mal: int=default_mal, mae: float=default_mae, aggmethod: str='median', filter: bool='both', outdir: str=None, no_shuffle: bool=False, nproc: int=-1):
        self.paired = bool(rev_format)
        assert fwd_format != rev_format, 'If "rev_format" is provided, it must differ from "fwd_format".\nOnly provide "fwd_format" for single end mode.'
        self.fwd_format = fwd_format
        self.rev_format = rev_format
        suffixes = tuple([fmt for fmt in [fwd_format, rev_format] if fmt])
        self.filepaths = sorted([f for f in filepaths if f.endswith(suffixes)])
        self.nreads = nreads
        self.mal = mal
        self.mae = mae
        assert aggmethod in ['mean', 'median'], '"aggmethod" must be either "median" (default) or "mean"'
        self.aggmethod = aggmethod
        self.aggfunc = np.median if aggmethod == 'median' else np.mean
        assert filter in ['any', 'both', 'no'], '"filter" must be either "any", "both", or "no"'
        self.filter = filter
        if outdir:
            self.outdir = os.path.abspath(outdir)
        else:
            self.outdir = None
        self.no_shuffle = no_shuffle
        self.nproc = nproc

        # All reads
        self.fwd_reads: List[Read] = []
        self.rev_reads : List[Read]= []
        self.samples: List[Sample] = None
        self._sample_mapping = None
        # Subset for minEER
        self.fwd_sub: List[Read] = None
        self.fwd_n_tries = None
        self.fwd_frac_passing = None
        self.rev_sub: List[Read] = None
        self.rev_n_tries = None
        self.rev_frac_passing = None
        self.all_sub = []
        # Global
        self.fwd_pos = None
        self.fwd_len = None
        self.rev_pos = None
        self.rev_len = None
        self.passing_readpairs: List[ReadPair] = None

    @functools.cached_property
    def files(self) -> List[File]:
        """All files in project"""
        fs = []
        for path in self.filepaths:
            direction = 'f' if path.endswith(self.fwd_format) else 'r'
            suffix = self.fwd_format if direction == 'f' else self.rev_format
            file = File(path, direction, self.mal, self.mae, suffix)
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
        Samples = [Sample(sample, files['f'], files['r'], self.filter) for sample, files in samples.items()]
        self.samples = Samples
        self._sample_mapping = samples

    @property
    def _reportInputs(self):
        """Report number of samples and reads"""
        if not self.samples:
            self.getReadsandSamples()
        # Basic inputs
        args = ['paired', 'fwd_format', 'rev_format', 'nreads', 'mal', 'mae', 'aggmethod', 'filter', 'outdir']
        nfiles = len(self.files)
        nsamples = len(self.samples)
        nreadpairs = len([rp for s in self.samples for rp in s.readpairs])
        pairs = {arg: getattr(self, arg) for arg in args}
        pairs.update({'Files': nfiles, 'Samples': nsamples, 'Read Pairs': nreadpairs})

        order = args + ['Files', 'Samples', 'Read Pairs']
        input_report = alignedSpacing(pairs, 50, order)

        # Sample-pairs
        sample_map_data = []
        for s, file_dict in self._sample_mapping.items():
            s_data = {'Sample': s, 'Forward': os.path.basename(file_dict['f'].filepath)}
            if self.paired:
                s_data.update({'Reverse': os.path.basename(file_dict['r'].filepath)})
            sample_map_data.append(s_data)
        sample_map_df = pd.DataFrame(sample_map_data).set_index('Sample')

        return input_report, sample_map_df
    
    def runMineer(self, reads: List[Read]):
        """Run mineer on a list of reads in parallel"""
        # Extract expected error profiles
        ees = [read.untrimmed.ee for read in reads]
        # Run through mineer
        trimposes = Parallel(self.nproc)(delayed(minEER)(ee, self.mal, self.mae) for ee in ees)
        # Update reads
        for read, trimpos in zip(reads, trimposes):
            read.trimpos_mineer = trimpos

    def subsampleReads(self, reads: List[Read]):
        """Subsample up to `nreads` that pass minEER"""
        if not self.no_shuffle: # ONLY EVER TRUE FOR TESTING
            # Shuffle reads
            random.shuffle(reads)

        # Keep track of passing
        read_count = 0
        passing_reads = []

        # Get passing reads until `nreads` is reached

        # Assume only 70% of reads in a chunk will pass
        multiplier = 1.3
        chunk_size = int(self.nreads * multiplier)
        start = 0
        while True:
            # Get positions to chunk
            end = start + chunk_size
            chunk = reads[start: end]
            # If an empty list, then no more reads to try
            if not chunk:
                break
            # Update number of reads tried
            read_count += len(chunk)
            # Run mineer on this chunk
            self.runMineer(chunk)
            # Collect passing
            curr_n_passing = 0
            for r in chunk:
                if r.pass_qc_mineer:
                    curr_n_passing += 1
                    passing_reads.append(r)
            # Calculate number of reads remaining
            remaining = self.nreads - len(passing_reads)
            # If we don't need anymore reads, break
            if remaining <= 0:
                break
            # Otherwise, get the next chunk with size(remaining)
            else:
                # Update multiplier by previous chunk's accuracy
                error = curr_n_passing / chunk_size
                multiplier = 1 + error
                chunk_size = int(remaining * multiplier)
                # Inefficient to have a chunk smaller than 1000
                if chunk_size < 1000:
                    chunk_size = 1000
                # Start off where we ended
                start = end

        return read_count, passing_reads
    
    def subsampleAll(self):
        """Subsample all reads up to `nreads` keeping those that pass minEER and recording fraction that pass"""
        self.fwd_n_tries, self.fwd_sub = self.subsampleReads(self.fwd_reads)
        self.fwd_frac_passing = len(self.fwd_sub) / self.fwd_n_tries
        self.all_sub.extend(self.fwd_sub)
        if self.paired:
            self.rev_n_tries, self.rev_sub = self.subsampleReads(self.rev_reads)
            self.rev_frac_passing = len(self.rev_sub) / self.rev_n_tries
            self.all_sub.extend(self.rev_sub)
    
    @property
    def _reportMineerStats(self):
        """Report fraction of passing reads from minEER run"""
        passing_stats = f'From subset run on minEER:\n\tForward reads passing: {len(self.fwd_sub)}/{self.fwd_n_tries} ({self.fwd_frac_passing * 100:.2f})%\n'
        if self.paired:
            passing_stats += f'\tReverse reads passing: {len(self.rev_sub)}/{self.rev_n_tries} ({self.rev_frac_passing * 100:.2f})%\n'

        """Report stats on trimpositions from minEER run"""
        fwd_poses = [{'trimstart': r.trimpos_mineer[0], 'trimend': r.trimpos_mineer[1], 'direction': 'f'} for r in self.fwd_sub if r.pass_qc_mineer]
        data = fwd_poses
        if self.paired:
            rev_poses = [{'trimstart': r.trimpos_mineer[0], 'trimend': r.trimpos_mineer[1], 'direction': 'r'} for r in self.rev_sub if r.pass_qc_mineer]
            data.extend(rev_poses)
        trimpos_df = pd.DataFrame(data)
        trimpos_stats = trimpos_df.groupby('direction').agg(['min', 'max', 'mean', 'median'])
        return passing_stats, trimpos_stats

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

    @property
    def _reportGlobalPos(self):
        """Report global truncation positions"""
        text = f"""Using {self.aggmethod} and {self.nreads} reads, truncation positions are:
        Forward: {self.fwd_pos}
        """
        if self.paired:
            text += f'Reverse: {self.rev_pos}'
        return text

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

    @property
    def _reportTruncStats(self):
        """Report passing reads after using global truncation positions"""
        report_dict = {
            'Forward reads': f'{len([r for r in self.fwd_reads if r.pass_qc])}/{len(self.fwd_reads)}',
            'Reverse reads': f'{len([r for r in self.rev_reads if r.pass_qc])}/{len(self.rev_reads)}',
            'Read pairs': f'{len(self.passing_readpairs)}/{len([rp for s in self.samples for rp in s.readpairs])}'
        }
        if self.paired:
            order = ['Forward reads', 'Reverse reads', 'Read pairs']
        else:
            order = ['Forward reads', 'Read pairs']

        report = 'Reads that pass quality control after truncation:\n'
        report += alignedSpacing(report_dict, order=order)
        return report

    def writeFile(self, file: File):
        """Write one file with truncated sequences. Use File to preserve original order"""
        # Generate output path
        outpath = os.path.join(self.outdir, f'{file.sample}_mineer{file.suffix}')
        # Create output directory if it doesn't exist
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        # Gather truncated reads where pairs pass qc
        with dnaio.FastqWriter(outpath) as fh:
            for r in file.reads:
                if r.bothpassing:
                    fh.write(r.trimmed.record)
            
    def writeFiles(self):
        """Write all truncated files"""
        for f in self.files:
            self.writeFile(f)

def alignedSpacing(pairs, width = 30, order=None):
    """Generate text with aligned spacing for each element of `pairs` of `width`"""
    text = ''
    if not order:
        order = pairs.keys()
    for k in order:
        v = str(pairs[k])
        spaces = ' ' * (width - len(k) - len(v))
        line = f'{k}:{spaces}{v}\n'
        text += line
    return text
