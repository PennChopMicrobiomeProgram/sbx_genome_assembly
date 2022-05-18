# sbx_WGS (Whole Genome Sequencing)

## Introduction

sbx_WGS is an extension for the sunbeam pipeline for de novo microbial genome assembly. This pipeline uses SPAdes for single genome assembly and CheckM and read map coverage of the assembled genome for quality assessment (sbx_SPARCQ.rules). In addition, this pipeline uses anvi'o for contamination evaluation and taxonomic assignment through single copy core genes, amongst other things (sbx_SCCG.rules). FASTQ files for Mycoplasma genitalium are provided for testing.

### Installation
1. Add the packages in sbx_WGS_env.yml to sunbeam environment.yml and install with your sunbeam's ./install.sh --update env

2. Add config.yml to sunbeam_config.yml
```
cat config.yml >> /path/to/sunbeam_config.yml
```
3. Recommended for cluster: include the memory specifications in cluster.json, especially for checkm_tree rule
4. For running larger assemblies on a cluster, export the TMPDIR variable so anvio knows what directory to use to set up its databases (see run_example.sh)

## Options for config.yml
threads (SPAdes, BWA, samtools) & cog_threads (anvio): # of threads to use for running these programs

checkm_yml (optional): YAML file containing a sample:rank and sample:taxon dictionary for CheckM parameters (see example);
rank is one of {life,domain,phylum,class,order,family,genus,species};
taxon is the taxonomic classification for the specified 'rank'

window_size (Read mapping): define the window size to do sliding window coverage

sampling (read mapping): define minimum length of the contig to slide over

metagenome: option for running this pipeline on metagenomic samples

## Analysis

The R markdown template to display the results from this pipeline and for analysis can be found in the Analysis_templates repository on Enterprise

## Issues with pipeline
2022-05-13:

rule concat_scgs that uses the `anvi-get-sequences-for-hmm-hits` command no longer works with MUSCLE version 5.1. This command still works with MUSCLE version 3.8 (this issue can be ignored if removing anvi'o from pipeline)

## TODO (see comments in sbx_SCCG.rules)

1. Remove anvi'o dependency
2. Add HMMER to identify SCCGs
3. Use number of SCCGs to evaluate contamination and completeness in genomes
4. Find method for taxonomic assignment for WGS
5. Remove rules for read mapping; modify sbx_mapping_withFilter to take assembled genomes from this pipeline

## Contributors
This extension was adopted from pipelines written by Scott Daniel, Jung-Jin Lee, Ceylan Tanes, and Louis Taylor. Spades rules were adopted from sbx_spades (https://github.com/sunbeam-labs/sbx_spades). Read mapping rules were adopted from the sunbeam pipeline (https://github.com/sunbeam-labs/sunbeam/tree/stable/rules/mapping) and sbx_mapping_withFilter (https://github.com/ctanes/sbx_mapping_withFilter). Anvio rules were adopted from anvio_pan (https://github.com/junglee0713/anvio_pan)