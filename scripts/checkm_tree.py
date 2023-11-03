import subprocess as sp
import sys
from pathlib import Path

with open(snakemake.log[0], "w") as log:
    args = ["checkm"]

    if (
        snakemake.params.checkm_yml
        and snakemake.wildcards.sample in snakemake.params.taxon_yml
    ):
        taxon = f"'{snakemake.params.taxon_yml[snakemake.wildcards.sample]}'"
        rank = str(snakemake.params.rank_yml[snakemake.wildcards.sample])
        args += ["taxonomy_wf", "-t", f"{snakemake.threads}", f"{rank}", f"{taxon}"]
    else:
        args += ["lineage_wf", "-t", f"{snakemake.threads}", "-x", "fasta"]

    args += [f"{snakemake.params.bins}", f"{snakemake.params.tree_output}"]

    try:
        checkm_output = sp.check_output(args)
    except sp.CalledProcessError as e:
        sys.stderr.write(e.output.decode())
        log.write(e.output.decode())
        sys.exit(e.returncode)
    sys.stdout.write(checkm_output.decode())
    log.write(checkm_output.decode())

    Path(snakemake.output).touch()
