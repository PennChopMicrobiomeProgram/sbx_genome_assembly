#!/usr/bin/env bash

#SBATCH --no-requeue
#SBATCH -n 1
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --output=slurm_%x_%j.out
#SBATCH -t 10-0

source ~/.bashrc.conda #needed to make "conda" command to work
conda activate sunbeam #needed for sunbeam

module load BioPerl

unset LD_LIBRARY_PATH #sunbeam / snakemake clears this path when it launched due to secrutiy issues

#technically, this should be set when conda activates the sunbeam environment
#but we set it manually just to be sure
export SUNBEAM_DIR="/home/tuv/sunbeam/sunbeam-stable"
#export TMPDIR variable for anvio steps
export TMPDIR="/prj/dir/tmp"
# -n after "all_decontam" does a dry run
# -n \

set -x
set -e

sunbeam run --configfile ./sunbeam_config.yml --use-conda all_WGS -j 50 -w90 -k --keep-going --ignore-incomplete --cluster-config ./cluster.json --notemp -p --nolock --verbose -c 'sbatch --no-requeue --export=ALL --mem={cluster.memcpu} -n {threads} -t 10-0 -J {cluster.name} --output=slurm_%x_%j.out'
