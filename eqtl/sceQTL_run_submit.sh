#!/bin/bash

###############################################################################
# snBD — Submit cell-type cis-eQTL jobs
#
# Purpose:
#   Submit one MatrixEQTL job per entry in data/cell_types_25.txt.
#
# Inputs:
#   Prepared inputs listed in sceQTL_run_each_cell_type.R; activated R.
#   data/cell_types_25.txt lists the cell types to analyze.
#
# Outputs:
#   out/sceQTL430/result/ and its logs directory.
#
# Method overview:
#   Read the selected cell types; create folders; submit the companion runner.
#
# Dependencies:
#   Bash, SLURM sbatch, Rscript and MatrixEQTL/genio dependencies.
#
# Usage:
#   From repository root: bash eqtl/sceQTL_run_submit.sh
#
# Before running:
#   Edit the cell-type list and scheduling parameters before submitting.
#   Prepare genotype, expression and covariate matrices with matching donor columns.
#
###############################################################################

#SBATCH --job-name test
#SBATCH -c 1                             
#SBATCH -p test
#SBATCH -t 10:00   
#SBATCH --mem=200M
#SBATCH -o ./logs/run_%x_%N_%j.out
#SBATCH -e ./logs/run_%x_%N_%j.err



# Pass the repository root to the analysis jobs.
# 1. Repository root and local environment -----------------------------
export SNBD_ROOT="${SNBD_ROOT:-$PWD}"

set -eu 


partition=shared
# 2. Resources and job parameters --------------------------------------
threads=1

mem=100000
walltime=2:00:00 

# module purge
# PATHTOSPACK=./software/spack
# source $PATHTOSPACK/share/spack/setup-env.sh
# spack load anaconda3

# source activate ./software/renv


dir_main="${SNBD_ROOT}"
# for cell_type in Ex-L23
# for cell_type in Endo


# 3. Submit jobs for the selected cell types ----------------------------
for cell_type in `awk 'NR >= 1 ' ${SNBD_ROOT}/data/cell_types_25.txt `
do


echo "${cell_type}"

sub_outpath="${dir_main}/out/sceQTL430/result"
log_path="${sub_outpath}/logs"


mkdir -p "${sub_outpath}"
mkdir -p "${log_path}"

mkdir -p "${dir_main}/out/sceQTL430/data"

	
sbatch --job-name "${cell_type}" --export=ALL,cell_type=${cell_type},sub_outpath=${sub_outpath},threads=${threads},mem=${mem} \
-e "${log_path}/${cell_type}.err" -o "${log_path}/${cell_type}.out" \
--cpus-per-task=${threads} --time=${walltime} --mem=${mem} --partition=${partition} \
"${SNBD_ROOT}/eqtl/sceQTL_run_each_cell_type.R"  --cell_type ${cell_type} --sub_outpath "${sub_outpath}";


sleep 1
echo -e "\n\n"

done


