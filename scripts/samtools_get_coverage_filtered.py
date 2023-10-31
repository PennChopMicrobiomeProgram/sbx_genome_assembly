from sunbeamlib import samtools

with open(snakemake.log[0], "w") as log:
    samtools.get_coverage_stats(
        snakemake.wildcards.sample,
        snakemake.input[0],
        snakemake.wildcards.sample,
        snakemake.output[0],
    )
