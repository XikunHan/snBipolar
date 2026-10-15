#!/bin/bash

###############################################################################
# snBD — Submit per-cell-type MR workers
#
# Purpose:
#   Submit one scMR worker per selected cell type.
#
# Inputs:
#   data/cell_types_25.txt; prepared eQTL and MR inputs.
#   Edit partition and resource settings below.
#
# Outputs:
#   out/sceQTL430/scMR/IV/clump/logs/ and worker outputs.
#
# Method overview:
#   Pass cell-type and path settings to scMR_run_each.sh and submit one job per cell type.
#
# Dependencies:
#   Bash, SLURM sbatch and the environment for scMR_run_each.sh.
#
# Usage:
#   From repository root: bash mendelian_randomization/scMR_run_submit.sh
#
###############################################################################

#SBATCH --job-name test
#SBATCH -c 1                             
#SBATCH -p test
#SBATCH -t 4:00:00   
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
walltime=10:00:00


# module load R/4.2.0


dir_main="${SNBD_ROOT}/out/sceQTL430"

# for cell_type in  Oli

# 3. Submit jobs for the selected cell types ----------------------------
for cell_type in `awk 'NR >= 1 ' ${SNBD_ROOT}/data/cell_types_25.txt `
do


echo "${cell_type}"

sub_outpath="${dir_main}/scMR/IV/clump"
log_path="${sub_outpath}/logs"


mkdir -p "${sub_outpath}"
mkdir -p "${log_path}"

	
sbatch --job-name "${cell_type}" --export=ALL,cell_type=${cell_type},sub_outpath=${sub_outpath},dir_main=${dir_main} \
-e "${log_path}/${cell_type}.err" -o "${log_path}/${cell_type}.out" \
--cpus-per-task=${threads} --time=${walltime} --mem=${mem} --partition=${partition} \
"${SNBD_ROOT}/mendelian_randomization/scMR_run_each.sh";


sleep 1
echo -e "\n\n"

done



