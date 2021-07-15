"""
minEER - A sequence trimming algorithm that seeks to minimize the expected error rate and maximize sequence length

AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

default_mal = 100
default_mae = 1e-3

class File:
    """Fastq file"""
    def __init__(self, filepath: str, mal: int=default_mal, mae: float=default_mae):
        self.filepath = filepath
        self.mal = mal
        self.mae = mae
    @property
    def reads(self):
        """Load all reads"""
        return [Read(r, self, self.mal, self.mae) for r in SeqIO.parse(self.filepath, 'fastq')]


class Record:
    """A single untrimmed read"""
    def __init__(self, seqrecord):
        self.record = seqrecord
        
    @property
    def phred(self):
        """Extract phred quality scores"""
        return self.record.letter_annotations['phred_quality']
    
    @property
    def ee(self):
        """Get expected error scores"""
        return phred2ee(self.phred)
    
class Read:
    """Untrimmed and Trimmed data on read"""
    def __init__(self, seqrecord: SeqRecord, file: File, mal: int=default_mal, mae: float=default_mae):
        self.untrimmed = Record(seqrecord)
        self.file = file
        self.mal = mal
        self.mae = mae
        self.pass_qc = None
        self.trimmed = None
    
    @property
    def trimpos(self):
        """Get trim positions using minEER"""
        return minEER(self.untrimmed.ee, self.mal, self.mae)
    
    def trim(self):
        """Trim read using minEER"""
        # Get trim positions using minEER
        trimstart, trimend = self.trimpos
        # Check if QC is passed 
        if (trimstart is None) & (trimend is None):
            self.pass_qc = False
        else:
            self.pass_qc = True
            # Trim sequence
            trimmed_seqrecord = self.untrimmed.record[trimstart: trimend + 1]
            self.trimmed = Record(trimmed_seqrecord)

def minEER(ee, mal=100, mae=1e-3):
    """
    Exhaustive expected error rate (EER) search for a sequence
    Finds the longest sequence with an EER<=thresh with length>=mal
    Inputs:
    | ee <array>: Vector of expected errors
    | mal <int>: Minimum acceptable length of sequence
    | mae <float>: Minimum acceptable EER threshold
    Outputs:
    | {start, end}_pos: Positions where subsequence starts and ends
    """
    # Ensure np.array()
    ee = np.array(ee)
    # Initialize matrix
    n_pos = len(ee)
    eer_space = np.empty((n_pos, n_pos))
    eer_space[:] = np.nan
    # For each starting position
    for start in range(n_pos - mal):
        # Get expected error profile from this start
        ee_start = ee[start:]
        # Calculate the expected errors starting at this position
        eers = ee_start.cumsum() / (np.arange(len(ee_start)) + 1)
        # Only take the expected errors starting at the minimum acceptable length
        eers_mal = eers[mal:]
        # Update the search space
        eer_space[start, (start + mal) :] = eers_mal

    """Find max length sequence below threshold"""
    # Filter out nans to avoid warning about trying to evaluate np.nan<thresh
    idxs = np.array(np.where(np.isfinite(eer_space)))
    threshed = eer_space[idxs[0], idxs[1]] <= mae
    # If no subsequences, return None
    if not threshed.any():
        start_pos, end_pos = None, None
    else:
        # Start and end indices that satisfy threshold
        start_ids, end_ids = idxs[:, threshed]
        # Find longest subsequence
        max_idx = np.argmax(end_ids - start_ids)
        start_pos, end_pos = start_ids[max_idx], end_ids[max_idx]
    return start_pos, end_pos 

def phred2ee(phred):
    """
    Convert a Phred score to an expected error probability
    ee = 10 ^ (-P/10)
    """
    phred = np.array(phred)
    ee = np.power(10, -phred / 10)
    return ee

# if __name__ == '__main__':
#     parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
#     parser.add_argument('-i', help='Input, type specified with -t', type=str, required=True)
#     parser.add_argument('-t', help='Input type:\n'
#                                    '\tfile: If -i is a path to a single file'
#                                    '\tdir: If -i is a path to a directory only containing files to trim'
#                                    '\tlist: If -i is a path to a list of paths of files to trim',
#                         choice=['file', 'dir', 'list'], type=str, required=True)
#     # TODO: allow per-file specificiation (could have input be CSV with filepath and format)
#     parser.add_argument('-f', help='File format of file(s) passed in -i', required=True)
#     parser.add_argument('-o', help='Output directory\n'
#                                    'Output files will be named as /output_dir/{original_name}_trimmed.{extension}')