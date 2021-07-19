"""
minEER - A sequence trimming algorithm that seeks to minimize the expected error rate and maximize sequence length

AUTHOR: MICHAEL SILVERSTEIN
EMAIL: michael.silverstein4@gmail.com
"""
import numpy as np
from numba import jit

@jit(nopython=True)
def eerspace(ee, mal):
    """
    Compute expected error rate space where each row is starting position and each column is an ending position

    For each starting position, calculate the expected error rate of every subsequence of the minimum acceptable length
    """
    # Initialize matrix
    n_pos = len(ee)
    eer_space = np.empty((n_pos, n_pos))
    eer_space[:] = np.nan
    # Compute the cumulative sum of the expected errors
    cumsum_ee = np.cumsum(ee)
    # Get normalizing vector to compute EER at each position
    norm = np.arange(n_pos) + 1
    # Initialize first start
    eer_space[0, mal:] = (cumsum_ee / norm)[mal:]

    # For each starting position
    for start in range(1, n_pos - mal):
        # Get cumulative sum from this start position
        cumsum_ee_start = cumsum_ee[start:]

        # Calculate the expected errors from this start position
        eers = (cumsum_ee_start - cumsum_ee[start - 1]) / norm[:-start]
        # Only take the expected errors starting at the minimum acceptable length
        eers_mal = eers[mal:]
        # Update the search space
        eer_space[start, (start + mal) :] = eers_mal
    return eer_space


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
    # Calculate expected error space
    eer_space = eerspace(ee, mal)

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
        # Add one to end pos to be inclusive when indexing
        end_pos += 1
    return start_pos, end_pos 