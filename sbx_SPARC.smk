# -*- mode: Snakemake -*-
#
# Rules for de novo assembly using SPAdes and post-assembly assessments

from sunbeamlib import samtools
import glob
import pysam
import re
import yaml


ruleorder: run_spades_paired > run_spades_unpaired


rule run_spades_paired:
    input:
        r1=str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
        r2=str(QC_FP / "decontam" / "{sample}_2.fastq.gz"),
    output:
        str(ASSEMBLY_FP / "spades_bins" / "{sample}" / "contigs.fasta"),
    benchmark:
        BENCHMARK_FP / "run_spades_paired_{sample}.tsv"
    log:
        LOG_FP / "run_spades_paired_{sample}.tsv",
    params:
        output_dir=ASSEMBLY_FP / "spades_bins" / "{sample}",
        log_fp=ASSEMBLY_FP / "spades_bins" / "{sample}" / "spades.log",
    threads: Cfg["sbx_WGS"]["threads"]
    conda:
        "sbx_WGS_env.yml"
    shell:
        """
        spades.py -1 {input.r1} -2 {input.r2} -o {params.output_dir} -t {threads} --cov-cutoff 5.0 --phred-offset 64 2>&1 | tee {log}; \
        cat {params.log_fp} >> {log}
        """


rule run_spades_unpaired:
    input:
        r1=str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
    output:
        str(ASSEMBLY_FP / "spades_bins" / "{sample}" / "{sample}_assembled_contigs.fna"),
    benchmark:
        BENCHMARK_FP / "run_spades_unpaired_{sample}.tsv"
    log:
        LOG_FP / "run_spades_unpaired_{sample}.tsv",
    params:
        outdir=str(ASSEMBLY_FP / "spades" / "{sample}"),
        mk_dir=str(ASSEMBLY_FP / "spades_bins" / "sample"),
        copy_from=str(ASSEMBLY_FP / "spades" / "{sample}" / "contigs.fasta"),
    threads: Cfg["sbx_WGS"]["threads"]
    conda:
        "sbx_WGS_env.yml"
    shell:
        """
        spades.py --s 1 {input.r1} -o {params.outdir} -t {threads} --cov-cutoff 5.0 --phred-offset 33 2>&1 | tee {log} && \
        mkdir -p {params.mk_dir} && \
        cp {params.copy_from} {output}
        """


checkm_yml = Cfg["sbx_WGS"].get("checkm_yml")
if checkm_yml:
    with open(checkm_yml, "r") as file:
        checkm_params = yaml.load(file, Loader=yaml.FullLoader)
        rank_yml = checkm_params["rank"]
        taxon_yml = checkm_params["taxon"]


rule checkm_tree:
    input:
        str(ASSEMBLY_FP / "spades_bins" / "{sample}" / "contigs.fasta"),
    output:
        str(ASSEMBLY_FP / "checkm_output" / "tree_output" / "{sample}" / "tree_done"),
    benchmark:
        BENCHMARK_FP / "checkm_tree_{sample}.tsv"
    log:
        LOG_FP / "checkm_tree_{sample}.log",
    params:
        bins=str(ASSEMBLY_FP / "spades_bins" / "{sample}"),
        tree_output=str(ASSEMBLY_FP / "checkm_output" / "tree_output" / "{sample}"),
        checkm_yml=checkm_yml,
        rank_yml=rank_yml if "rank_yml" in locals() else None,
        taxon_yml=taxon_yml if "taxon_yml" in locals() else None,
    threads: Cfg["sbx_WGS"]["threads"]
    conda:
        "sbx_WGS_env.yml"
    script:
        "scripts/checkm_tree.py"


rule checkm_summary:
    input:
        str(ASSEMBLY_FP / "checkm_output" / "tree_output" / "{sample}" / "tree_done"),
    output:
        str(
            ASSEMBLY_FP
            / "checkm_output"
            / "summary"
            / "{sample}"
            / "extended_summary.tsv"
        ),
    benchmark:
        BENCHMARK_FP / "checkm_summary_{sample}.tsv"
    log:
        LOG_FP / "checkm_summary_{sample}.log",
    params:
        tree_output=str(ASSEMBLY_FP / "checkm_output" / "tree_output" / "{sample}"),
        checkm_yml=checkm_yml,
        taxon_yml=taxon_yml if "taxon_yml" in locals() else None,
    conda:
        "sbx_WGS_env.yml"
    script:
        "scripts/checkm_summary.py"


def _checkm_summary_tsvs(w):
    pattern = str(
        ASSEMBLY_FP / "checkm_output" / "summary" / "{sample}" / "extended_summary.tsv"
    )
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule checkm_summarize:
    input:
        _checkm_summary_tsvs,
    output:
        str(ASSEMBLY_FP / "checkm_output" / "all_extended_summary.tsv"),
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"


rule index_assembled_genomes:
    input:
        str(ASSEMBLY_FP / "anvio" / "{sample}" / "{sample}_reformatted_contigs.fa"),  #this file comes from anvio
    output:
        str(
            ASSEMBLY_FP
            / "read_mapping"
            / "{sample}"
            / "bwa"
            / "{sample}_reformatted_contigs.fa.amb"
        ),
    benchmark:
        BENCHMARK_FP / "index_assembled_genomes_{sample}.tsv"
    log:
        LOG_FP / "index_assembled_genomes_{sample}.log",
    params:
        bwa_dir=str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "bwa"),
        bwa_sample=str(
            ASSEMBLY_FP
            / "read_mapping"
            / "{sample}"
            / "bwa"
            / "{sample}_reformatted_contigs.fa"
        ),
    conda:
        "sbx_WGS_env.yml"
    shell:
        """
        mkdir -p {params.bwa_dir} && \
        cp {input} {params.bwa_dir} && \
        cd {params.bwa_dir} && \
        bwa index {params.bwa_sample} 2>&1 | tee {log}
        """


rule align_2_genome:
    input:
        reads=expand(str(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz"), rp=Pairs),
        index=str(
            ASSEMBLY_FP
            / "read_mapping"
            / "{sample}"
            / "bwa"
            / "{sample}_reformatted_contigs.fa.amb"
        ),
    output:
        temp(
            str(
                ASSEMBLY_FP
                / "read_mapping"
                / "{sample}"
                / "bwa"
                / "intermediates"
                / "{sample}.sam"
            )
        ),
    benchmark:
        BENCHMARK_FP / "align_2_genome_{sample}.tsv"
    log:
        LOG_FP / "align_2_genome_{sample}.log",
    params:
        str(
            ASSEMBLY_FP
            / "read_mapping"
            / "{sample}"
            / "bwa"
            / "{sample}_reformatted_contigs.fa"
        ),
    threads: Cfg["sbx_WGS"]["threads"]
    conda:
        "sbx_WGS_env.yml"
    shell:
        """
        bwa mem -M -t {threads} \
        {params} \
        {input.reads} -o {output} \
        2>&1 | tee {log}
        """


rule assembly_samtools_convert:
    input:
        str(
            ASSEMBLY_FP
            / "read_mapping"
            / "{sample}"
            / "bwa"
            / "intermediates"
            / "{sample}.sam"
        ),
    output:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "bwa" / "{sample}.bam"),
    benchmark:
        BENCHMARK_FP / "assembly_samtools_convert_{sample}.tsv"
    log:
        view_log=LOG_FP / "assembly_samtools_convert_view_{sample}.log",
        sort_log=LOG_FP / "assembly_samtools_convert_sort_{sample}.log",
    threads: Cfg["sbx_WGS"]["threads"]
    conda:
        "sbx_WGS_env.yml"
    shell:
        """
        samtools view -@ {threads} -b {input} 2>&1 | tee {log.view_log} | \
        samtools sort -@ {threads} -o {output} 2>&1 | tee {log.sort_log}
        """


rule index_samtools:
    input:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "bwa" / "{sample}.bam"),
    output:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "bwa" / "{sample}.bam.bai"),
    benchmark:
        BENCHMARK_FP / "index_samtools_{sample}.tsv"
    log:
        LOG_FP / "index_samtools_{sample}.log",
    conda:
        "sbx_WGS_env.yml"
    shell:
        "samtools index {input} {output} 2>&1 | tee {log}"


rule samtools_get_coverage_filtered:
    input:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "bwa" / "{sample}.bam"),
    output:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "genome_coverage_{sample}.csv"),
    benchmark:
        BENCHMARK_FP / "samtools_get_coverage_filtered_{sample}.tsv"
    log:
        LOG_FP / "samtools_get_coverage_filtered_{sample}.log",
    conda:
        "sbx_WGS_env.yml"
    script:
        "scripts/samtools_get_coverage_filtered.py"


def _sorted_csvs(w):
    pattern = str(
        ASSEMBLY_FP / "read_mapping" / "{sample}" / "genome_coverage_{sample}.csv"
    )
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule summarize_assembly_coverage:
    input:
        _sorted_csvs,
    output:
        str(ASSEMBLY_FP / "read_mapping" / "all_coverage.csv"),
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"


rule samtools_summarize_num_mapped_reads:
    input:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "bwa" / "{sample}.bam"),
    output:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "numReads_{sample}.csv"),
    benchmark:
        BENCHMARK_FP / "samtools_summarize_num_mapped_reads_{sample}.tsv"
    log:
        LOG_FP / "samtools_summarize_num_mapped_reads_{sample}.log",
    conda:
        "sbx_WGS_env.yml"
    shell:
        """
        samtools idxstats {input} | tee {log} | (sed 's/^/{wildcards.sample}\t/') > {output}
        """


def _numReads(w):
    pattern = str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "numReads_{sample}.csv")
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule samtools_summarize_numReads:
    input:
        _numReads,
    output:
        str(ASSEMBLY_FP / "read_mapping" / "all_numReads.csv"),
    shell:
        "(cat {input}) > {output}"


def sliding_window_coverage(genome, bamfile, sample, output_fp, N, sampling):
    print(genome)
    print(bamfile)
    print(sample)
    print(output_fp)
    print(N)
    print(sampling)

    output_rows = []
    args = ["samtools", "depth", "-aa", bamfile]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)
    # Organize into a list of depths for each segment, streaming in text
    reader = csv.reader(p.stdout, delimiter="\t")
    data = {}
    for row in reader:
        if not data.get(row[0]):
            data[row[0]] = []
        data[row[0]].append(int(row[2]))

    fields = ["Genome", "Segment", "Sample", "Location", "Average"]
    with open(output_fp, "w") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        for segment in data.keys():
            if len(data[segment]) > sampling:
                moving_avg = numpy.convolve(
                    data[segment], numpy.ones((N,)) / N, mode="full"
                )
                for i, x in enumerate(moving_avg):
                    if i % sampling == 0:
                        writer.writerow([genome, segment, sample, i, x])


rule samtools_get_sliding_coverage:
    input:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "bwa" / "{sample}.bam"),
    output:
        str(ASSEMBLY_FP / "read_mapping" / "{sample}" / "sliding_coverage_{sample}.csv"),
    benchmark:
        BENCHMARK_FP / "samtools_get_sliding_coverage_{sample}.tsv"
    log:
        LOG_FP / "samtools_get_sliding_coverage_{sample}.log",
    params:
        window_size=Cfg["sbx_WGS"]["window_size"],
        sampling=Cfg["sbx_WGS"]["sampling"],
        sliding_window_coverage=sliding_window_coverage,
    conda:
        "sbx_WGS_env.yml"
    script:
        "scripts/samtools_get_sliding_coverage.py"


def _sliding_coverage_csvs(w):
    pattern = str(
        ASSEMBLY_FP / "read_mapping" / "{sample}" / "sliding_coverage_{sample}.csv"
    )
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule samtools_summarize_sliding_coverage:
    input:
        _sliding_coverage_csvs,
    output:
        str(ASSEMBLY_FP / "read_mapping" / "all_sliding_coverage.csv"),
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"
