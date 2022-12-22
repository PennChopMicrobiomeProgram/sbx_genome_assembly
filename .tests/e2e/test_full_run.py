import csv
import os
import pytest
import subprocess as sp


@pytest.fixture(scope="session")
def dir(pytestconfig):
    return pytestconfig.getoption("dir")

@pytest.fixture
def setup(dir):
    temp_dir = dir

    reads_fp = os.path.abspath(".tests/data/reads/")

    project_dir = os.path.join(temp_dir, "project/")

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    yield temp_dir, project_dir


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    # Run the test job
    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_WGS",
            "--directory",
            temp_dir,
        ]
    )

    output_fp = os.path.join(project_dir, "sunbeam_output/")

    all_SCCG_fp = os.path.join(output_fp, "assembly/hmmer/all_SCCG_hits.tsv")

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield all_SCCG_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    all_SCCG_fp, benchmarks_fp = run_sunbeam

    with open(all_SCCG_fp) as f:
        assert len(f.readlines()) > 140, f"Wasn't able to find at least 70 hits"


def test_benchmarks(run_sunbeam):
    all_SCCG_fp, benchmarks_fp = run_sunbeam

    filename = os.listdir(benchmarks_fp)[0]
    with open(os.path.join(benchmarks_fp, filename)) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            assert (
                float(r["cpu_time"]) < 0.5
            ), f"cpu_time for {r['rule']} is higher than 0.5: {r['cpu_time']}"
