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

        self.samples_fp = os.path.join(self.temp_dir, "samples.csv")
        self.samples_content = """
ASM2732v1,{r0r1},{r0r2}\n
Mixed,{r1r1},{r1r2}\n
""".format(
                r0r1 = os.path.join(self.data_fp, "ASM2732v1_fixed_R1.fastq.gz"),
                r0r2 = os.path.join(self.data_fp, "ASM2732v1_fixed_R2.fastq.gz"),
                r1r1 = os.path.join(self.data_fp, "Mixed_fixed_R1.fastq.gz"),
                r1r2 = os.path.join(self.data_fp, "Mixed_fixed_R2.fastq.gz"))
        with open(self.samples_fp, "w") as f:
            f.write(self.samples_content)
        
        self.config_fp = os.path.join(self.temp_dir, "sunbeam_config.yml")
        self.config_content = """
            all:\n
              root: {root}\n
              output_fp: sunbeam_output\n
              samplelist_fp: samples.csv\n
              paired_end: true\n
              download_reads: false\n
              version: 2.1.1+dev0.g43432e1.d20220208\n
            qc:\n
              suffix: qc\n
              seq_id_ending: ''\n
              threads: 4\n
              java_heapsize: 512M\n
              leading: 3\n
              trailing: 3\n
              slidingwindow: [4, 15]\n
              minlen: 36\n
              adapter_fp: {adapt}\n
              adapter_template: $CONDA_PREFIX/share/trimmomatic/adapters/NexteraPE-PE.fa\n
              fwd_adapters: [GTTTCCCAGTCACGATC, GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC]\n
              rev_adapters: [GTTTCCCAGTCACGATC, GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC]\n
              kz_threshold: 0.55\n
              pct_id: 0.5\n
              frac: 0.6\n
              host_fp: ''\n
            classify:\n
              suffix: classify\n
              threads: 4\n
              kraken_db_fp: ''\n
            assembly:\n
              suffix: assembly\n
              min_length: 300\n
              threads: 4\n
            annotation:\n
              suffix: annotation\n
              min_contig_len: 500\n
              circular_kmin: 10\n
              circular_kmax: 1000\n
              circular_min_len: 3500\n
            blast:\n
              threads: 4\n
            blastdbs:\n
              root_fp: ''\n
            mapping:\n
              suffix: mapping\n
              genomes_fp: ''\n
              samtools_opts: ''\n
              threads: 4\n
            download:\n
              suffix: download\n
              threads: 4\n
              sbx_coassembly: ''\n
              threads: 4\n
              group_file: ''\n
            sbx_WGS:\n
              threads: 4\n
              cog_db:\n
              cog_threads: 4\n
              metagenome: false\n
              checkm_yml:\n
              window_size: 5000\n
              sampling: 5000\n
            """.format(root = os.path.join(self.temp_dir),
                adapt = os.path.join(find_conda(), "envs/sunbeam/share/trimmomatic/adapters/NexteraPE-PE.fa"))
        with open(self.config_fp, "w") as f:
            f.write(self.config_content)

        self.output_fp = os.path.join(self.temp_dir, "sunbeam_output")

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