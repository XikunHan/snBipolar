#!/bin/bash

###############################################################################
# snBD — Submit cell-type SMR jobs
#
# Purpose:
#   Submit one SMR worker per selected cell type.
#
# Inputs:
#   data/cell_types_25.txt; prepared BESD, LD-reference and GWAS files.
#   Edit scheduler resources and binary locations as needed.
#
# Outputs:
#   out/sceQTL430/SMR/logs/ and worker outputs in SMR_out/.
#
# Method overview:
#   Pass cell-type and path settings to SMR_run_each.sh and submit one job per cell type.
#
# Dependencies:
#   Bash, SLURM sbatch, SMR and the worker environment.
#
# Usage:
#   From repository root: bash smr/SMR_run_submit.sh
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
mem=20000
walltime=1:00:00


# module load R/4.2.0

dir_main="${SNBD_ROOT}"
# for cell_type in In-Rosehip_CHST9
# 3. Submit jobs for the selected cell types ----------------------------
for cell_type in `awk 'NR >= 1 ' ${SNBD_ROOT}/data/cell_types_25.txt `
do


echo "${cell_type}"

sub_outpath="${dir_main}/out/sceQTL430/SMR"
log_path="${sub_outpath}/logs"

mkdir -p "${sub_outpath}"
mkdir -p "${log_path}"

sbatch --job-name "${cell_type}" --export=ALL,cell_type=${cell_type},sub_outpath=${sub_outpath},threads=${threads},mem=${mem} \
-e "${log_path}/${cell_type}_multi.err" -o "${log_path}/${cell_type}_multi.out" \
--cpus-per-task=${threads} --time=${walltime} --mem=${mem} --partition=${partition} \
"${SNBD_ROOT}/smr/SMR_run_each.sh";


sleep 1
echo -e "\n\n"

done


