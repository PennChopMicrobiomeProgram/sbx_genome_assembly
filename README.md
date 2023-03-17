<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_genome_assembly

[![Tests](https://github.com/PennChopMicrobiomeProgram/sbx_genome_assembly/actions/workflows/main.yml/badge.svg)](https://github.com/PennChopMicrobiomeProgram/sbx_genome_assembly/actions/workflows/main.yml)
[![Super-Linter](https://github.com/PennChopMicrobiomeProgram/sbx_genome_assembly/actions/workflows/linter.yml/badge.svg)](https://github.com/PennChopMicrobiomeProgram/sbx_genome_assembly/actions/workflows/linter.yml)

## Introduction

sbx_genome_assembly is an extension for the sunbeam pipeline for de novo microbial genome assembly. This pipeline uses SPAdes for single genome assembly and CheckM and read map coverage of the assembled genome for quality assessment (sbx_SPARC rules). In addition, this pipeline uses hmmer to identify SCCGs (sbx_SCCG.rules). FASTQ files for Mycoplasma genitalium are provided for testing.

### Installation

```
sunbeam extend https://github.com/PennChopMicrobiomeProgram/sbx_genome_assembly.git
```

Recommended for cluster: include the memory specifications in cluster.json, especially for checkm_tree rule

## Running

Run with sunbeam on the target `all_genome_assembly`,

```
sunbeam run --configfile /path/to/sunbeam_config.yml all_genome_assembly
```

### Options for config.yml
threads (SPAdes, BWA, samtools): # of threads to use for running these programs

checkm_yml (optional): YAML file containing a sample:rank and sample:taxon dictionary for CheckM parameters (see example);
rank is one of {life,domain,phylum,class,order,family,genus,species};
taxon is the taxonomic classification for the specified 'rank'

window_size (Read mapping): define the window size to do sliding window coverage

sampling (read mapping): define minimum length of the contig to slide over

metagenome: option for running this pipeline on metagenomic samples

## Analysis

The R markdown template to display the results from this pipeline and for analysis can be found in the Analysis_templates repository on Enterprise

## Legacy Installation

```
git clone https://github.com/PennChopMicrobiomeProgram/sbx_genome_assembly.git extensions/sbx_genome_assembly
cd extensions/sbx_genome_assembly
cat config.yml >> /path/to/sunbeam_config.yml
```

## Issues with pipeline

Please post any issues with this extension [here](https://github.com/PennChopMicrobiomeProgram/sbx_genome_assembly/issues).

## TODO (see comments in sbx_SCCG.rules)

1. Use number of SCCGs to evaluate contamination and completeness in genomes
2. Find method for taxonomic assignment for genome assembly
3. Remove rules for read mapping; modify sbx_mapping_withFilter to take assembled genomes from this pipeline

## Contributors
This extension was adopted from pipelines written by Scott Daniel, Jung-Jin Lee, Ceylan Tanes, and Louis Taylor. Spades rules were adopted from sbx_spades (https://github.com/sunbeam-labs/sbx_spades). Read mapping rules were adopted from the sunbeam pipeline (https://github.com/sunbeam-labs/sunbeam/tree/stable/rules/mapping) and sbx_mapping_withFilter (https://github.com/ctanes/sbx_mapping_withFilter). Anvio rules were adopted from anvio_pan (https://github.com/junglee0713/anvio_pan)
