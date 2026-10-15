#!/bin/bash

###############################################################################
# snBD — Submit spatial-registration jobs
#
# Purpose:
#   Submit one spatialLIBD registration job for McLean and one for MtSinai.
#
# Inputs:
#   Prepared combined SCE and an activated R/spatialLIBD environment.
#   Edit SLURM partition/resources below.
#
# Outputs:
#   out/sp_registration/<cohort>/ and logs for each cohort.
#
# Method overview:
#   Create output/log folders and call sbatch_sp_registration.R.
#
# Dependencies:
#   Bash, SLURM sbatch, Rscript and spatial-registration R packages.
#
# Usage:
#   From repository root: bash spatial_registration/sbatch_sp_registration.sh
#
###############################################################################

#SBATCH --job-name module
#SBATCH -c 1                             
#SBATCH -p serial_requeue
#SBATCH -t 10:00   
#SBATCH --mem=10M
#SBATCH -o ./log/run_%x_%N_%j.out
#SBATCH -e ./log/run_%x_%N_%j.err



# Pass the repository root to the analysis jobs.
# 1. Repository root and local environment -----------------------------
export SNBD_ROOT="${SNBD_ROOT:-$PWD}"

set -eu 


# module purge
# PATHTOSPACK=./software/spack
# source $PATHTOSPACK/share/spack/setup-env.sh
# spack load anaconda3

# source activate ./software/renv4.3




# 2. Resources and job parameters --------------------------------------
threads=1
mem=100000 
walltime=4:00:00
partition=shared



dir_main="${SNBD_ROOT}/out/sp_registration"



for cohort in  McLean  MtSinai
do


sub_outpath="${dir_main}/${cohort}"
log_path="${sub_outpath}/logs"


mkdir -p "${sub_outpath}"
mkdir -p "${log_path}"


sbatch --job-name "sp_${cohort}" --export=ALL \
-e "${log_path}/${cohort}.err" -o "${log_path}/${cohort}.out" \
--cpus-per-task=${threads} --time=${walltime} --mem=${mem} --partition=${partition} \
"${SNBD_ROOT}/spatial_registration/sbatch_sp_registration.R"  --sub_outpath "${sub_outpath}" --cohort ${cohort};


sleep 1
echo -e "\n\n"

done

