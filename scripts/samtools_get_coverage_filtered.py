from sunbeamlib import samtools

samtools.get_coverage_stats(snakemake.wildcards.sample, snakemake.input[0], snakemake.wildcards.sample, snakemake.output[0])