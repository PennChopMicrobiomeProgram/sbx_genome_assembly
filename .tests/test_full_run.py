import os
import sys

import subprocess as sp
import shutil
import unittest
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

def find_conda() -> str:
    for fp in sys.path:
        if "miniconda" in fp or "anaconda" in fp:
            conda_fp = ""
            for s in fp.split("/"):
                conda_fp = conda_fp + s + "/"
                if "miniconda" in s or "anaconda" in s:
                    print("Found conda installation at " + conda_fp)
                    return conda_fp
    print("WARNING: Couldn't find path to NexteraPE adapter")

class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.data_fp = os.path.join(self.temp_dir, "data/")
        shutil.copytree("test_samples/", self.data_fp)
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
            "summarize_hmmscan",
            "--directory",
            self.temp_dir,
        ])

        # Check output
        self.assertTrue(os.path.exists(self.all_SCCG_fp))
        with open(self.all_SCCG_fp) as f:
            self.assertGreater(len(f.readlines()), 140) # Found at least 70 SCCGs in each
