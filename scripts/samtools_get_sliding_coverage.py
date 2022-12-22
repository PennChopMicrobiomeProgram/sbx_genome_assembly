import sys

with open(snakemake.log[0], "w") as l:
    sys.stderr = sys.stdout = l
    snakemake.params.sliding_window_coverage(
        snakemake.wildcards.sample,
        snakemake.input[0],
        snakemake.wildcards.sample,
        snakemake.output[0],
        snakemake.params.window_size,
        snakemake.params.sampling,
    )
