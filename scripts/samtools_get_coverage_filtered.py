import sys
from sunbeamlib import samtools

with open(snakemake.log[0], "w") as l:
    sys.stderr = sys.stdout = l
    samtools.get_coverage_stats(
        snakemake.wildcards.sample,
        snakemake.input[0],
        snakemake.wildcards.sample,
        snakemake.output[0],
    )
