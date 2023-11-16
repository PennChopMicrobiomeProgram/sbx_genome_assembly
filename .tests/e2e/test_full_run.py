import csv
import os
import pytest
import shutil
import subprocess as sp
import tempfile
from pathlib import Path


@pytest.fixture
def setup():
    temp_dir = Path(tempfile.mkdtemp())

    reads_fp = Path(".tests/data/reads/").resolve()
    hosts_fp = Path(".tests/data/hosts/").resolve()

    project_dir = temp_dir / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"qc: {{host_fp: {hosts_fp}}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    # Run the test job.
    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "test_genome_assembly",
            "--directory",
            temp_dir,
        ]
    )

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    contigs_fp = output_fp / "assembly/spades_bins/GCA_0004638950/contigs.fasta"
    all_SCCG_fp = output_fp / "assembly/hmmer/all_SCCG_hits.tsv"

    assert os.path.exists(contigs_fp)
    assert os.stat(contigs_fp).st_size > 0
    assert os.path.exists(all_SCCG_fp)
    with open(all_SCCG_fp) as f:
        assert len(f.readlines()) > 140, f"Wasn't able to find at least 70 hits"
