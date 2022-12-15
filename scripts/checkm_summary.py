import os

if (
    snakemake.params.checkm_yml
    and snakemake.wildcards.sample in snakemake.params.taxon_yml
):
    taxon = str(snakemake.params.taxon_yml[snakemake.wildcards.sample])
    os.system(
        """
    checkm qa --out_format 2 --tab_table --file {snakemake.output} "{snakemake.params.tree_output}/{taxon}.ms" {snakemake.params.tree_output}
    """
    )
else:
    os.system(
        """
    checkm qa --out_format 2 --tab_table --file {snakemake.output} "{snakemake.params.tree_output}/lineage.ms" {snakemake.params.tree_output}
    """
    )
