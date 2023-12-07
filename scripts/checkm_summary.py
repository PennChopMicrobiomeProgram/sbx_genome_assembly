import subprocess as sp
import sys

with open(snakemake.log[0], "w") as log:
    if (
        snakemake.params.checkm_yml
        and snakemake.wildcards.sample in snakemake.params.taxon_yml
    ):
        taxon = str(snakemake.params.taxon_yml[snakemake.wildcards.sample])

    try:
        checkm_output = sp.check_output(
            [
                "checkm",
                "qa",
                "--out_format",
                "2",
                "--tab_table",
                "--file",
                f"{snakemake.output}",
                f"{snakemake.input.lineage}",
                f"{snakemake.params.tree_output}",
            ]
        )
    except sp.CalledProcessError as e:
        log.write(e.output.decode())
        sys.exit(e.returncode)
    log.write(checkm_output.decode())
