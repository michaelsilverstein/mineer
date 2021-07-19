"""Test pipeline"""

from Bio import SeqIO
from unittest import TestCase
from mineer.pipeline import truncPipeline
import os, subprocess, shutil, random

class TestPipeline(TestCase):
    """Test whole pipeline on two samples"""
    @classmethod
    def setUpClass(self):
        self.tempdir = 'test/test_files/fastqs'
        # Download fastqs
        subprocess.call(f'fastq-dump SRR9660307 SRR9660321 --split-files -O {self.tempdir}'.split(' '))
        self.filepaths = [f'{os.path.join(self.tempdir, f)}' for f in os.listdir(self.tempdir)]

        self.fwd_format = '_1.fastq'
        self.rev_format = '_2.fastq'

        self.project = truncPipeline(self.filepaths, self.fwd_format, self.rev_format, outdir=self.tempdir, random_seed=123)
    
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.tempdir)
    
    def test_n_tries(self):
        # Forward
        self.assertEqual(self.project.fwd_n_tries, 5128)
        # Reverse
        self.assertEqual(self.project.rev_n_tries, 7223)
    
    def test_trimpos(self):
        self.assertEqual(self.project.fwd_pos.tolist(), [0, 301])
        self.assertEqual(self.project.rev_pos.tolist(), [0, 164])

    def test_passing_readpairs(self):
        self.assertEqual(self.project.passing_readpairs, 1872)

    def test_reads_same_len(self):
        self.assertTrue(all([r.trimmed.length == self.project.fwd_len for r in self.project.fwd_reads if r.pass_qc]))
        self.assertTrue(all([r.trimmed.length == self.project.rev_len for r in self.project.rev_reads if r.pass_qc]))
        
    def test_writing_readpairs(self):
        for sample in self.project.samples:
            # Saved to file
            r1_fastqs = len(list(SeqIO.parse(f'{self.tempdir}/{sample.name}_mineer_1.fastq', 'fastq')))
            r2_fastqs = len(list(SeqIO.parse(f'{self.tempdir}/{sample.name}_mineer_2.fastq', 'fastq')))
            sample_passing_readpairs = len([rp for rp in sample.readpairs if rp.bothpassing])
            self.assertTrue(r1_fastqs == r2_fastqs == sample_passing_readpairs)


