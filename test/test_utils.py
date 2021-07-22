"""
Test mineer utils
"""
from mineer.utils import Project, ReadPair, Read
from Bio.SeqRecord import SeqRecord
from unittest import TestCase
from copy import copy

class TestProjectPaired(TestCase):
    def setUp(self):
        samples = ['a', 'b']
        suffixes = ['_1.fastq', '_2.fastq']
        filepaths = [f'{sample}{suffix}' for sample in samples for suffix in suffixes]
        self.kwargs = dict(
            filepaths = filepaths,
            fwd_format = suffixes[0],
            rev_format = suffixes[1],
            nreads = 10,
            mal = 20,
            mae = 1,
            aggmethod = 'mean',
            outdir = 'test'
        )
        self.project = Project(**self.kwargs)

    def test_project_ingest(self):
        for arg, val in self.kwargs.items():
            project_val = getattr(self.project, arg)
            self.assertEqual(val, project_val)
    
    def test_ignore_bad_files(self):
        bad_paths = ['bad', 'paths']
        expected = copy(self.kwargs['filepaths'])
        self.kwargs['filepaths'].extend(bad_paths)
        self.assertListEqual(self.project.filepaths, expected)

class TestProjectSingle(TestCase):
    def setUp(self):
        filepaths = ['a.fastq', 'b.fastq']
        self.kwargs = dict(
        filepaths = filepaths,
        fwd_format = '.fastq',
        nreads = 10,
        mal = 20,
        mae = 1,
        aggmethod = 'mean',
        outdir = 'test'
        )
        self.project = Project(**self.kwargs)

    def test_project_ingest(self):
        for arg, val in self.kwargs.items():
            project_val = getattr(self.project, arg)
            self.assertEqual(val, project_val)
    
    def test_not_paired(self):
        self.assertFalse(self.project.paired)

    def test_only_forward(self):
        bad_paths = ['a_2.fastq', 'b_2.fastq']
        expected = copy(self.kwargs['filepaths'])
        self.kwargs['filepaths'].extend(bad_paths)
        self.assertListEqual(self.project.filepaths, expected)

class FakeRead:
    # Hack the system
    def __init__(self, pass_l, pass_e):
        self.trimmed = True
        self.pass_l = pass_l
        self.pass_e = pass_e
        self.pass_qc = self.pass_e & self.pass_l

class TestReadPair(TestCase):
    def test_filter_any(self):
        # Bothpassing case
        r1 = FakeRead(True, True)
        r2 = FakeRead(True, True)
        rp = ReadPair(r1, r2, filter='any')
        self.assertTrue(rp.bothpassing)

        # Failing
        r2 = FakeRead(True, False)
        rp = ReadPair(r1, r2, 'any')
        self.assertFalse(rp.bothpassing)
    
    def test_filter_both(self):
        # Bothpassing case
        r1 = FakeRead(True, True)
        r2 = FakeRead(True, False)
        rp = ReadPair(r1, r2, filter='both')
        self.assertTrue(rp.bothpassing)

        # Failing
        r1 = FakeRead(True, False)
        rp = ReadPair(r1, r2, 'both')
        self.assertFalse(rp.bothpassing)

    def test_filter_no(self):
        # Bothpassing case
        r1 = FakeRead(True, False)
        r2 = FakeRead(True, False)
        rp = ReadPair(r1, r2, filter='no')
        self.assertTrue(rp.bothpassing)

        # Failing
        r1 = FakeRead(False, True)
        r2 = FakeRead(False, False)
        rp = ReadPair(r1, r2, 'no')
        self.assertFalse(rp.bothpassing)     