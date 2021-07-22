"""Test pipeline"""

from mineer.download_test_files import download
from Bio import SeqIO
from unittest import TestCase
from mineer.pipeline import truncPipeline, mineer_cli
import os, shutil
from copy import copy

class TestPipeline(TestCase):
    """Test whole pipeline on two samples"""
    @classmethod
    def setUpClass(self):
        self.tempdir = 'test/test_files/fastqs'
        self.out = 'test/test_files/out'
        self.filter_out = 'test/test_files/filter'
        self.cli_out = 'test/test_files/cli_out'
        self.viz_out = 'test/test_files/test_viz'
        self.single_out = 'test/test_files/single'
        # Download fastqs
        download(self.tempdir)
        self.kwargs = dict(
            outdir = self.out,
            filepaths = [f'{os.path.join(self.tempdir, f)}' for f in os.listdir(self.tempdir)],
            fwd_format = '_1.fastq',
            rev_format = '_2.fastq',
            nreads = 1000,
            no_shuffle = True
        )
    
        self.project = truncPipeline(**self.kwargs, viz_outdir=self.viz_out)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.out)
        shutil.rmtree(self.cli_out)
        shutil.rmtree(self.viz_out)
        shutil.rmtree(self.single_out)
    
    def test_n_tries(self):
        # Forward
        self.assertEqual(self.project.fwd_n_tries, 1030)
        # Reverse
        self.assertEqual(self.project.rev_n_tries, 1643)
    
    def test_trimpos(self):
        self.assertEqual(self.project.fwd_pos.tolist(), [0, 301])
        self.assertEqual(self.project.rev_pos.tolist(), [0, 166])

    def test_passing_readpairs(self):
        self.assertEqual(len(self.project.passing_readpairs), 4975)

    def test_reads_same_len(self):
        self.assertTrue(all([r.trimmed.length == self.project.fwd_len for r in self.project.fwd_reads if r.pass_qc]))
        self.assertTrue(all([r.trimmed.length == self.project.rev_len for r in self.project.rev_reads if r.pass_qc]))
        
    def test_writing_readpairs(self):
        """Test that same number of reads are in r1 and r2 as are passing"""
        for sample in self.project.samples:
            outdir = self.kwargs['outdir']
            # Saved to file
            r1_fastqs = len(list(SeqIO.parse(f'{outdir}/{sample.name}_mineer_1.fastq', 'fastq')))
            r2_fastqs = len(list(SeqIO.parse(f'{outdir}/{sample.name}_mineer_2.fastq', 'fastq')))
            sample_passing_readpairs = len([rp for rp in sample.readpairs if rp.bothpassing])
            self.assertTrue(r1_fastqs == r2_fastqs == sample_passing_readpairs)

    def test_cli(self):
        cmd = f'-i {self.tempdir} -f _1.fastq -r _2.fastq -n 1000 -o {self.cli_out} --test'
        cli_project = mineer_cli(cmd.split())
        p_rps = len([rp for s in self.project.samples for rp in s.readpairs])
        c_rps = len([rp for s in cli_project.samples for rp in s.readpairs])
        self.assertEqual(p_rps, c_rps)
        self.assertEqual(len(cli_project.passing_readpairs), len(self.project.passing_readpairs))
        self.assertListEqual(cli_project.fwd_pos.tolist(), self.project.fwd_pos.tolist())
        self.assertListEqual(cli_project.rev_pos.tolist(), cli_project.rev_pos.tolist())
        
        # Check files are there
        outfiles = sorted(os.listdir(self.cli_out))
        self.assertEqual(outfiles, ['SRR9660307_mineer_1.fastq', 'SRR9660307_mineer_2.fastq', 'SRR9660321_mineer_1.fastq', 'SRR9660321_mineer_2.fastq'])

    def test_viz(self):
        self.assertEqual(sorted(os.listdir(self.viz_out)), ['phred_profiles.png', 'trunc_dist.png'])
    
    def test_filter_no(self):
        kwargs = copy(self.kwargs)
        kwargs['outdir'] = self.filter_out
        project = truncPipeline(**kwargs, filter='no')
        self.assertEqual(len(project.passing_readpairs), 6741)

    def test_filter_any(self):
        kwargs = copy(self.kwargs)
        kwargs['outdir'] = self.filter_out
        project = truncPipeline(**kwargs, filter='any')
        self.assertEqual(len(project.passing_readpairs), 1812)
    
    def test_single_end(self):
        kwargs = copy(self.kwargs)
        kwargs.pop('rev_format')
        kwargs['outdir'] = self.single_out
        project = truncPipeline(**kwargs)
        # No 'rev' args
        rev_arg_vals = [v for k, v in project.__dict__.items() if 'rev' in k]
        self.assertFalse(any(rev_arg_vals))

        outfiles = sorted(os.listdir(self.single_out))
        self.assertEqual(outfiles, ['SRR9660307_mineer_1.fastq', 'SRR9660321_mineer_1.fastq'])
