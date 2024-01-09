#!/bin/bash
#
#SBATCH --job-name=test_METARes
#SBATCH --mail-user=pierre.lemee@anses.fr
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=160:00:00

. /local/env/envconda.sh
. /local/env/envsnakemake-7.28.3.sh

snakemake -s METARes.snakefile -k --use-conda -j 32 --latency-wait 60 --conda-frontend mamba