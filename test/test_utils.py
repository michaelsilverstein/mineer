"""
Test mineer utils
"""
from mineer.utils import Project
from unittest import TestCase
from copy import copy

class TestProject(TestCase):
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
    
    def test_project_ingest(self):
        project = Project(**self.kwargs)
        for arg, val in self.kwargs.items():
            project_val = getattr(project, arg)
            self.assertEqual(val, project_val)
    
    def test_ignore_bad_files(self):
        bad_paths = ['bad', 'paths']
        expected = copy(self.kwargs['filepaths'])
        self.kwargs['filepaths'].extend(bad_paths)
        project = Project(**self.kwargs)
        self.assertListEqual(project.filepaths, expected)

    def test_not_paired(self):
        self.kwargs.pop('rev_format')
        project = Project(**self.kwargs)
        self.assertFalse(project.paired)

