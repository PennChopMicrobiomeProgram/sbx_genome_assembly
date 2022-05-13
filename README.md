# sbx_WGS (Whole Genome Sequencing)

## Introduction

sbx_WGS is an extension for the sunbeam pipeline for de novo microbial genome assembly and quality assessment. This pipeline uses SPAdes for single genome assembly and CheckM, QUAST, and Read map coverage of the assembled genome for quality assessment (sbx_SPARCQ.rules). In addition, it uses anvi'o for contamination assessment and taxonomic assignment with single copy genes, amongst other things (sbx_SCG.rules). An R markdown template has been provided that will nicely display the results generated from this pipeline. Illumina reads are provided for testing.

### Installation
1. Add packages in sbx_SPARCQ_env.yml to sunbeam environment.yml and install with your sunbeam's ./install.sh --update env
2. To install QUAST in sunbeam environment (if you cannot install it in step #1), use your environment's pip to install quast
```
/path/to/miniconda3/envs/sunbeam/bin/pip install quast
(If installed this way, may need to change quast to quast.py in sbx_SPARCQ.rules)
```
3. Add config.yml to sunbeam_config.yml
```
cat config.yml >> /path/to/sunbeam_config.yml
```
4. Recommended for cluster: add the memory specifications in cluster.json, especially for checkm_tree rule
5. For running larger assemblies on a cluster, export the TMPDIR variable so anvio knows what directory to use to set up its databases (see run_example.sh)

## Options for config.yml
threads (SPAdes, BWA, samtools) & cog_threads (anvio): # of threads to use for running programs

checkm_yml (optional): YAML file containing a sample:rank and sample:taxon dictionary for CheckM parameters (see example);
rank is one of {life,domain,phylum,class,order,family,genus,species};
taxon is the taxonomic classification for the specified 'rank'

taxid_yml (optional): YAML file containing a sample:TaxonID dictionary for reference genomes to be downloaded by ncbi-genome-download for comparison in QUAST

window_size (Read mapping): define the window size to do sliding window coverage

sampling (read mapping): define minimum length of the contig to slide over

metagenome: option for running this pipeline on metagenomic samples

## Analysis

Templates for analysis can be found in the Analysis_templates repository

## Issues with pipeline
2022-05-13:

rule concat_scgs that uses the `anvi-get-sequences-for-hmm-hits` command no longer works with MUSCLE version 5.1. This command still works with MUSCLE version 3.8

rule scg_setup which uses `anvi-setup-scg-databases` fails to download the require database. This may be due to GTDB updating its database on April 2022 and changing its path to the db

## Contributors
This extension was adapted from pipelines written by Scott Daniel, Jung-Jin Lee, Ceylan Tanes, and Louis Taylor. Spades rules were adapted from sbx_spades (https://github.com/sunbeam-labs/sbx_spades). Read mapping rules were adapted from the sunbeam pipeline (https://github.com/sunbeam-labs/sunbeam/tree/stable/rules/mapping) and sbx_mapping_withFilter (https://github.com/ctanes/sbx_mapping_withFilter). Anvio rules were adapted from anvio_pan (https://github.com/junglee0713/anvio_pan)