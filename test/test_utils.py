"""
Test mineer utils
"""
from mineer.utils import Project
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

    def test_project_ingest(self):
        project = Project(**self.kwargs)
        for arg, val in self.kwargs.items():
            project_val = getattr(project, arg)
            self.assertEqual(val, project_val)
    
    def test_not_paired(self):
        project = Project(**self.kwargs)
        self.assertFalse(project.paired)

    def test_only_forward(self):
        bad_paths = ['a_2.fastq', 'b_2.fastq']
        expected = copy(self.kwargs['filepaths'])
        self.kwargs['filepaths'].extend(bad_paths)
        self.assertListEqual(self.project.filepaths, expected)
