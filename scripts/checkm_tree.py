import os
import sys

with open(snakemake.log[0], "w") as l:
    sys.stderr = sys.stdout = l
    if (
        snakemake.params.checkm_yml
        and snakemake.wildcards.sample in snakemake.params.taxon_yml
    ):
        taxon = str(snakemake.params.taxon_yml[snakemake.wildcards.sample])
        rank = str(snakemake.params.rank_yml[snakemake.wildcards.sample])
        os.system(
            f"""
        checkm taxonomy_wf -t {snakemake.threads} {rank} '{taxon}' {snakemake.params.bins} {snakemake.params.tree_output} && \
        touch {snakemake.output}
        """
        )
    else:
        os.system(
            f"""
        checkm lineage_wf -t {snakemake.threads} {snakemake.params.bins} {snakemake.params.tree_output} && \
        touch {snakemake.output}
        """
        )
