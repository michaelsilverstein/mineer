"""Test minEER algorithm"""

from unittest import TestCase
import numpy as np
import numba
numba.config.DISABLE_JIT = True
from mineer.mineer import eerspace, minEER
from mineer.utils import Record
from Bio import SeqIO

class TestToy(TestCase):
    def setUp(self):
        self.ee = [1, 2, 3]
        self.mal = 0
        self.mae = np.inf

    def test_eerspace(self):
        # Compute eerspace
        space = eerspace(self.ee, self.mal)
        # Expected eerspace
        expected = [[1., 1.5, 2.],
                    [np.nan, 2., 2.5],
                    [np.nan, np.nan, 3.]]
        comp = np.array_equal(space, expected, equal_nan=True)
        self.assertTrue(comp)
        
    def test_mineer(self):
        # Get trim positions
        trimpos = minEER(self.ee, self.mal, self.mae)
        self.assertEqual(trimpos, (0, 3))

class TestFile(TestCase):
    def setUp(self):
        record = Record(SeqIO.read('test/test_files/test_read', 'fastq'))
        self.ee = record.ee
    
    def test_truncpos(self):
        truncpos = minEER(self.ee)
        expected = (80, 206)
        self.assertEqual(truncpos, expected)
        