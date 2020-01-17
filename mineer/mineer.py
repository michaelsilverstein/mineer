"""
minEER - A sequence trimming algorithm that seeks to minimize the expected error rate and maximize sequence length

AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
import warnings
import numpy as np
from argparse import ArgumentParser, RawTextHelpFormatter
from Bio.SeqRecord import SeqRecord

class TrimObject:
    """Class to manage sequence information"""
    def __init__(self, Seq, trimstart, trimend):
        """
        Initialize SeqObject
        """
        self.untrimmed = Seq
        self.trimstart = None
        self.trimend = None
        self.trimmed = None


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

def trim(Seq, mal=100, mae=1e-3):
    """
    Trims sequence object `Seq` according to parameters `mal` and `mae`
    Inputs:
    | Seq <Bio.SeqRecord>: Sequence object
    | mal <int>: Minimum acceptable length
    | mae <float>: Minimum acceptable expected error rate
    Outputs:
    | Trimmed <Bio.Seq>: Trimmed sequence object
    """
    # Check for SeqRecord
    if not isinstance(Seq, SeqRecord):
        raise TypeError('"Seq" must be a BioPython SeqRecord object.\n'
                        'Read the docs on how to load a SeqRecord here: https://biopython.org/wiki/SeqIO')
    # Check to see that sequence contains PHRED score
    annot = Seq.letter_annotations
    if 'phred_qualty' in annot:
        phred = annot['phred_quality']
        # Convert to expected error
        ee = phred2ee(phred)
        # Implement minEER to get trimming indices
        trimstart, trimmed = minEER(ee, mal, mae)
        # Trim
        trimmed = Seq[trimstart: trimmed]
    else:
        warnings.warn('Sequence "%s" does not contain a phred score - '
                      'it cannot be trimmed and will be returned as is.' % Seq.id)
        trimmed = Seq
    return trimmed

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', help='Input, type specified with -t', type=str, required=True)
    parser.add_argument('-t', help='Input type:\n'
                                   '\tfile: If -i is a path to a single file'
                                   '\tdir: If -i is a path to a directory only containing files to trim'
                                   '\tlist: If -i is a path to a list of paths of files to trim',
                        choice=['file', 'dir', 'list'], type=str, required=True)
    # TODO: allow per-file specificiation (could have input be CSV with filepath and format)
    parser.add_argument('-f', help='File format of file(s) passed in -i', required=True)
    parser.add_argument('-o', help='Output directory\n'
                                   'Output files will be named as /output_dir/{original_name}_trimmed.{extension}')