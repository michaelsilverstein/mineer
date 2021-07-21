"""Test pipeline"""

from mineer.utils import Project
from mineer.download_test_files import download
from Bio import SeqIO
from unittest import TestCase
from mineer.pipeline import truncPipeline, mineer_cli
import os, shutil

class TestPipelinePaired(TestCase):
    """Test whole pipeline on two samples"""
    @classmethod
    def setUpClass(self):
        self.tempdir = 'test/test_files/fastqs'
        self.outdir1 = 'test/test_files/test1'
        self.outdir2 = 'test/test_files/test2'
        self.viz_out = 'test/test_files/test_viz'
        # Download fastqs
        download(self.tempdir)
        self.filepaths = [f'{os.path.join(self.tempdir, f)}' for f in os.listdir(self.tempdir)]

        self.fwd_format = '_1.fastq'
        self.rev_format = '_2.fastq'
        self.nreads = 5000

        self.project = truncPipeline(self.filepaths, self.fwd_format, self.rev_format, nreads=self.nreads, outdir=self.outdir1, random_seed=123)
    
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.outdir1)
        shutil.rmtree(self.outdir2)
        shutil.rmtree(self.viz_out)
    
    # def test_n_tries(self):
    #     # Forward
    #     self.assertEqual(self.project.fwd_n_tries, 5128)
    #     # Reverse
    #     self.assertEqual(self.project.rev_n_tries, 7223)
    
    # def test_trimpos(self):
    #     self.assertEqual(self.project.fwd_pos.tolist(), [0, 301])
    #     self.assertEqual(self.project.rev_pos.tolist(), [0, 164])

    def test_passing_readpairs(self):
        self.assertEqual(len(self.project.passing_readpairs), 1872)

    def test_reads_same_len(self):
        self.assertTrue(all([r.trimmed.length == self.project.fwd_len for r in self.project.fwd_reads if r.pass_qc]))
        self.assertTrue(all([r.trimmed.length == self.project.rev_len for r in self.project.rev_reads if r.pass_qc]))
        
    def test_writing_readpairs(self):
        """Test that same number of reads are in r1 and r2 as are passing"""
        for sample in self.project.samples:
            # Saved to file
            r1_fastqs = len(list(SeqIO.parse(f'{self.outdir1}/{sample.name}_mineer_1.fastq', 'fastq')))
            r2_fastqs = len(list(SeqIO.parse(f'{self.outdir1}/{sample.name}_mineer_2.fastq', 'fastq')))
            sample_passing_readpairs = len([rp for rp in sample.readpairs if rp.bothpassing])
            self.assertTrue(r1_fastqs == r2_fastqs == sample_passing_readpairs)

    def test_cli(self):
        cmd = f'-i {self.tempdir} -f _1.fastq -r _2.fastq -n 1000 -o {self.outdir2}'
        mineer_cli(cmd.split())
        # Check files are there
        outfiles = sorted(os.listdir(self.outdir2))
        self.assertEqual(outfiles, ['SRR9660307_mineer_1.fastq', 'SRR9660307_mineer_2.fastq', 'SRR9660321_mineer_1.fastq', 'SRR9660321_mineer_2.fastq'])

    def test_viz(self):
        project = truncPipeline(self.filepaths, self.fwd_format, self.rev_format, nreads=1000, outdir=self.outdir1, random_seed=123, viz_outdir=self.viz_out)
        self.assertEqual(os.listdir(self.viz_out), ['phred_profiles.png', 'trunc_dist.png'])

class TestPipelineSingle(TestCase):
    """Uses files downloaded from Paired"""
    @classmethod
    def setUpClass(self):
        self.tempdir = self.tempdir = 'test/test_files/fastqs'
        self.filepaths = [f'{os.path.join(self.tempdir, f)}' for f in os.listdir(self.tempdir)]
        self.project = Project(self.filepaths, '_1.fastq')


    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.tempdir)

    def test_forward_args_only(self):
        """All reverse arguments should be empty"""
        rev_arg_vals = [v for k, v in self.project.__dict__.items() if 'rev' in k]
        self.assertFalse(any(rev_arg_vals))