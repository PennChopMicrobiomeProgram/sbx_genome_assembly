# -*- mode: Snakemake -*-
#
# Rules for de novo assembly using SPAdes and post-assembly assessments

import sys

try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


localrules:
    all_genome_assembly,
    test_genome_assembly,


rule all_genome_assembly:
    input:
        # annotated ORFs from assembled genomes
        # genome quality summary from checkm
        # hmmscan hits on SCCGs
        [
            expand(
                ANNOTATION_FP / "prokka" / "{sample}" / "{sample}.faa",
                sample=Samples.keys(),
            ),
            ASSEMBLY_FP / "checkm_output" / "all_extended_summary.tsv",
            ASSEMBLY_FP / "hmmer" / "all_SCCG_hits.tsv",
        ],


rule test_genome_assembly:
    input:
        # skip checkm, it takes too much memory for GHActions
        [
            expand(
                ANNOTATION_FP / "prokka" / "{sample}" / "{sample}.faa",
                sample=Samples.keys(),
            ),
            ASSEMBLY_FP / "hmmer" / "all_SCCG_hits.tsv",
        ],


def get_input(wildcards):
    if Cfg["sbx_genome_assembly"]["metagenome"]:
        return ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa"
    else:
        return ASSEMBLY_FP / "spades_bins" / "{sample}" / "contigs.fasta"


rule reformat_fasta:
    input:
        get_input,
    output:
        str(ASSEMBLY_FP / "hmmer" / "{sample}" / "{sample}_reformatted_contigs.fa"),
    benchmark:
        BENCHMARK_FP / "reformat_fasta_{sample}.tsv"
    log:
        LOG_FP / "reformat_fasta_{sample}.log",
    params:
        len=200,
    script:
        "scripts/reformat_fasta.py"


rule prokka:
    input:
        str(ASSEMBLY_FP / "hmmer" / "{sample}" / "{sample}_reformatted_contigs.fa"),
    output:
        str(ANNOTATION_FP / "prokka" / "{sample}" / "{sample}.faa"),
    benchmark:
        BENCHMARK_FP / "prokka_{sample}.tsv"
    log:
        LOG_FP / "prokka_{sample}.log",
    params:
        outdir=str(ANNOTATION_FP / "prokka" / "{sample}"),
    conda:
        "sbx_SCCG_prokka_env.yml"
    shell:
        """
        prokka --compliant --centre CHOP --outdir {params.outdir} --locustag {wildcards.sample} --prefix {wildcards.sample} --force {input} 2>&1 | tee {log}
        """


def get_genome_assembly_path() -> str:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_genome_assembly":
            return fp
    sys.exit(
        "Filepath for genome_assembly not found, are you sure it's installed under extensions/sbx_genome_assembly?"
    )


rule hmmpress:
    input:
        os.path.join(get_genome_assembly_path(), "genes.hmm"),
    output:
        expand(
            os.path.join(get_genome_assembly_path(), "genes.hmm.h3{suffix}"),
            suffix={"f", "i", "m", "p"},
        ),
    benchmark:
        BENCHMARK_FP / "hmmpress.tsv"
    log:
        LOG_FP / "hmmpress.log",
    conda:
        "sbx_SCCG_env.yml"
    shell:
        "hmmpress {input} 2>&1 | tee {log}"


rule hmmscan:
    input:
        faa=str(ANNOTATION_FP / "prokka" / "{sample}" / "{sample}.faa"),
        aux=expand(
            os.path.join(get_genome_assembly_path(), "genes.hmm.h3{suffix}"),
            suffix={"f", "i", "m", "p"},
        ),
    output:
        str(ASSEMBLY_FP / "hmmer" / "{sample}" / "{sample}_SCCG_hits.tsv"),
    benchmark:
        BENCHMARK_FP / "hmmscan_{sample}.log"
    log:
        LOG_FP / "hmmscan_{sample}.log",
    params:
        hmm=os.path.join(get_genome_assembly_path(), "genes.hmm"),
    conda:
        "sbx_SCCG_env.yml"
    shell:
        "hmmscan -o {log} --tblout {output} {params.hmm} {input.faa}"


rule summarize_hmmscan:
    input:
        expand(
            str(ASSEMBLY_FP / "hmmer" / "{sample}" / "{sample}_SCCG_hits.tsv"),
            sample=Samples.keys(),
        ),
    output:
        str(ASSEMBLY_FP / "hmmer" / "all_SCCG_hits.tsv"),
    shell:
        "cat {input} > {output}"
