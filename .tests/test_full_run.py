import os
import sys

import subprocess as sp
import shutil
import unittest
import tempfile

class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.data_fp = os.path.join(self.temp_dir, "data/")
        shutil.copytree(".tests/test_samples/", self.data_fp)
        os.remove(os.path.join(self.data_fp, "commands2generateReads.txt"))

        self.project_dir = os.path.join(self.temp_dir, "project/")

        sp.check_output([
            "sunbeam",
            "init",
            "--data_fp",
            self.data_fp,
            self.project_dir
        ])
        
        self.config_fp = os.path.join(self.project_dir, "sunbeam_config.yml")

        self.output_fp = os.path.join(self.project_dir, "sunbeam_output")

        self.all_SCCG_fp = os.path.join(self.output_fp, "assembly/hmmer/all_SCCG_hits.tsv")
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_full_run(self):
        # Run the test job.
        sp.check_output([
            "sunbeam",
            "run",
             "--configfile", 
            self.config_fp,
            "--use-conda",
            "-c4",
            "test_WGS",
            "--directory",
            self.temp_dir,
        ])

        # Check output
        self.assertTrue(os.path.exists(self.all_SCCG_fp))
        with open(self.all_SCCG_fp) as f:
            self.assertGreater(len(f.readlines()), 140) # Found at least 70 SCCGs in each
