#!/bin/bash

###############################################################################
# snBD — Submit CellChat jobs by cohort
#
# Purpose:
#   Submit the McLean and Sinai CellChat analyses to SLURM.
#
# Inputs:
#   Combined SCE listed in cellchat_run_each.R; activated R environment.
#   Edit partition, memory, walltime and output directory below.
#
# Outputs:
#   out/cellchat/<cohort>/major/ and cohort-specific logs.
#
# Method overview:
#   Submit one cellchat_run_each.R job per cohort using the resources below.
#
# Dependencies:
#   Bash, sbatch and the R environment required by the companion runner.
#
# Usage:
#   From repository root: bash cell_cell_communication/cellchat_run_submit.sh
#
#   Invoke this launcher with bash.
#
###############################################################################

#SBATCH --job-name cellchat
#SBATCH -c 1                             
#SBATCH -t 10:00   
#SBATCH --mem=10M
#SBATCH -o ./log/run_%x_%N_%j.out
#SBATCH -e ./log/run_%x_%N_%j.err



# set -eu 


# module purge
# PATHTOSPACK=./software/spack
# source $PATHTOSPACK/share/spack/setup-env.sh
# spack load anaconda3

# source activate ./software/renv
# module load gcc



# Pass the repository root to the analysis jobs.
# 1. Repository root and local environment -----------------------------
export SNBD_ROOT="${SNBD_ROOT:-$PWD}"

# 2. Resources and job parameters --------------------------------------
threads=1
mem=100000 
walltime=4:00:00
partition=shared

dir_main="${SNBD_ROOT}/out/cellchat"


# 3. Submit one job per cohort -----------------------------------------
for cell_type in major
do

for cohort in McLean Sinai
do


sub_outpath="${dir_main}/${cohort}/${cell_type}"
log_path="${sub_outpath}/logs"


mkdir -p "${sub_outpath}"
mkdir -p "${log_path}"


sbatch --job-name "cellchat${cell_type}_${cohort}" --export=ALL \
-e "${log_path}/${cell_type}_${cohort}.err" -o "${log_path}/${cell_type}_${cohort}.out" \
--cpus-per-task=${threads} --time=${walltime} --mem=${mem} --partition=${partition} \
"${SNBD_ROOT}/cell_cell_communication/cellchat_run_each.R"  --cell_type ${cell_type}  --sub_outpath "${sub_outpath}" --cohort ${cohort};


sleep 1
echo -e "\n\n"

done
done

