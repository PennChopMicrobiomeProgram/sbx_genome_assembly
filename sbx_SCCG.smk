# -*- mode: Snakemake -*-
#
# Rules for de novo assembly using SPAdes and post-assembly assessments

from sunbeamlib import samtools
from Bio import pairwise2, SeqIO
import glob
import pysam
import re
import yaml
import os
import sys


localrules:
    all_WGS,
    test_WGS,


rule all_WGS:
    input:
        # annotated ORFs from assembled genomes
        # genome quality summary from checkm
        # hmmscan hits on SCCGs
        [
            expand(
                str(ANNOTATION_FP / "prokka" / "{sample}" / "{sample}.faa"),
                sample=Samples.keys(),
            ),
            str(ASSEMBLY_FP / "checkm_output" / "all_extended_summary.tsv"),
            str(ASSEMBLY_FP / "hmmer" / "all_SCCG_hits.tsv"),
        ],


rule test_WGS:
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
    if Cfg["sbx_WGS"]["metagenome"]:
        return str(ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa")
    else:
        return str(
            ASSEMBLY_FP / "spades_bins" / "{sample}" / "contigs.fasta"
        )


rule reformat_fasta:
    input:
        get_input,
    output:
        str(ASSEMBLY_FP / "hmmer" / "{sample}" / "{sample}_reformatted_contigs.fa"),
    params:
        len=200,
    script:
        "scripts/reformat_fasta.py"


rule prokka:
    input:
        str(ASSEMBLY_FP / "hmmer" / "{sample}" / "{sample}_reformatted_contigs.fa"),
    output:
        str(ANNOTATION_FP / "prokka" / "{sample}" / "{sample}.faa"),
    params:
        outdir=str(ANNOTATION_FP / "prokka" / "{sample}"),
    conda:
        "sbx_SCCG_env.yml"
    shell:
        """
        prokka --compliant --centre CHOP --outdir {params.outdir} --locustag {wildcards.sample} --prefix {wildcards.sample} --force {input} 
        """


def get_WGS_path() -> str:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_WGS":
            return fp
    raise Error(
        "Filepath for WGS not found, are you sure it's installed under extensions/sbx_WGS?"
    )


rule hmmpress:
    input:
        os.path.join(get_WGS_path(), "genes.hmm"),
    output:
        expand(
            os.path.join(get_WGS_path(), "genes.hmm.h3{suffix}"),
            suffix={"f", "i", "m", "p"},
        ),
    conda:
        "sbx_SCCG_env.yml"
    shell:
        "hmmpress {input}"


rule hmmscan:
    input:
        faa=str(ANNOTATION_FP / "prokka" / "{sample}" / "{sample}.faa"),
        aux=expand(
            os.path.join(get_WGS_path(), "genes.hmm.h3{suffix}"),
            suffix={"f", "i", "m", "p"},
        ),
    output:
        str(ASSEMBLY_FP / "hmmer" / "{sample}" / "{sample}_SCCG_hits.tsv"),
    log:
        "logs/{sample}_hmmscan.log",
    params:
        hmm=os.path.join(get_WGS_path(), "genes.hmm"),
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
