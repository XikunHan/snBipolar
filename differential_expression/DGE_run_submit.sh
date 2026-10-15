#!/bin/bash


################################################################################
# snBD — SLURM launcher for per–cell-type × cohort DE jobs
#
# Purpose:
#   Submit edgeR, NEBULA or MAST-RE jobs by cell type and cohort.
#   Each job calls DGE_run_each_cell_type.R.
#
# Dataset:
#   snRNA-seq PFC (BD vs CON) across two cohorts (McLean, MtSinai).
#   Split per cell type × cohort using the companion script: create_cell_type_file.R.
#
# Analysis design:
#   - Run differential expression separately in each cohort;
#     aggregate results via fixed-effect meta-analysis downstream.
#   - Include RUVr factors (K=3) and per-cell QC covariates as specified
#     in the R runner.
#
# Inputs (from prior step):
#   <data_root>/cell_type_data_McLean/snBD_sce_format_<celltype>.rds
#   <data_root>/cell_type_data_Sinai/snBD_sce_format_<celltype>.rds
#
# Outputs (per job):
#   <dir_main>/<DGE>_<Cohort>_<method>/

#
# Resource profile (defaults; may be tuned per method below):
#   edgeR:                   1 CPU, 10G RAM, 1h
#   NEBULA:                  1 CPU, 100G RAM, 10h
#   MAST-lmmBlock:            1 CPU, 50G RAM, 4h; block 1 is selected below
#
# Software/Modules:
#   R >= 4.2.0 (module load R/4.2.0)
#
# Logging:
#   Job names and log filenames include the cohort, cell type and gene block.
#
# Usage:
#   Edit the parameter section (trait/method/partition/paths/method list),
#   then run:
#     bash differential_expression/DGE_run_submit.sh
#
# Date: 2023-3
################################################################################



# module purge
# PATHTOSPACK=./software/spack
# source $PATHTOSPACK/share/spack/setup-env.sh
# spack load anaconda3
# source activate ./software/renv



# Pass the repository root to the analysis jobs.
# 1. Repository root and local environment -----------------------------
export SNBD_ROOT="${SNBD_ROOT:-$PWD}"

trait="BD"        # BD vs CON
method="RUV"     # RUV method label for outputs
partition="shared"


# Methods to run (comment/uncomment as needed)
# declare -a DGE_METHODS=(edger nebula)
# declare -a DGE_METHODS=(edger)
declare -a DGE_METHODS=(MASTlmmBlock)  # MAST with donor random effects

# declare -a COHORTS=(McLean Sinai)
declare -a COHORTS=(McLean)

# 25 cell types 
declare -a CELL_TYPES=(
  Ast Endo Ex-L23 Ex-L2 Ex-L3 Ex-L45_LRRK1 Ex-L45_MET Ex-L4_MYLK
  Ex-L56_CC_NTNG2 Ex-L56 Ex-L5b_HTR2C Ex-L6b_SEMA3D Ex-L6b_SEMA3E Ex-L6_CC_SEMA3A
  In-PV_Basket In-PV_Chandelier In-Reelin In-Rosehip_CHST9 In-Rosehip_TRPC3 In-SST
  In-VIP Mic Oli OPC Pericytes
)


# Test run: Endo only. Remove this assignment to use all 25 cell types.
declare -a CELL_TYPES=(
  Endo 
)



dir_main="${SNBD_ROOT}/out/DGE/test"
runner="${SNBD_ROOT}/differential_expression/DGE_run_each_cell_type.R"



for DGE in "${DGE_METHODS[@]}"; do
  for cohort in "${COHORTS[@]}"; do

    # defaults per method
    # 2. Resources and job parameters --------------------------------------
    threads=1
    walltime="1:00:00"
    mem_gb="10G"
    if [[ "$DGE" == "nebula" ]]; then
      walltime="10:00:00"
      mem_gb="100G"
    elif [[ "$DGE" == "MASTlmmBlock" ]]; then
      walltime="4:00:00"
      mem_gb="50G"
    fi

    sub_outpath="${dir_main}/${DGE}_${cohort}_${method}"
    log_path="${sub_outpath}/logs"
    mkdir -p "${sub_outpath}" "${log_path}"

    # 3. Submit jobs for the selected cell types ----------------------------
    for cell_type in "${CELL_TYPES[@]}"; do
      ct_safe="$(echo "$cell_type" | tr ' /' '__')"

      if [[ "$DGE" == "MASTlmmBlock" ]]; then
        # for block in $(seq 1 68); do
        for block in 1; do
          jobname="b${block}_${DGE}_${cohort}_${method}_${ct_safe}_${trait}"
          err="${log_path}/${trait}__${ct_safe}__b${block}.err"
          out="${log_path}/${trait}__${ct_safe}__b${block}.out"

         sbatch \
            --job-name="${jobname}" \
            --cpus-per-task="${threads}" \
            --time="${walltime}" \
            --mem="${mem_gb}" \
            --partition="${partition}" \
            -e "${err}" -o "${out}" \
            --export=ALL \
            --wrap "Rscript --vanilla \"$runner\" --trait ${trait} --DGE ${DGE} --cohort ${cohort} --method ${method} --cell_type ${cell_type} --sub_outpath \"${sub_outpath}\" --block ${block}"
          sleep 0.2
        done

      else
        # Non-block methods use block=NA.
        block="NA"
        jobname="${DGE}_${cohort}_${method}_${ct_safe}_${trait}"
        err="${log_path}/${trait}__${ct_safe}.err"
        out="${log_path}/${trait}__${ct_safe}.out"

        sbatch \
          --job-name="${jobname}" \
          --cpus-per-task="${threads}" \
          --time="${walltime}" \
          --mem="${mem_gb}" \
          --partition="${partition}" \
          -e "${err}" -o "${out}" \
          --export=ALL \
         --wrap "Rscript --vanilla \"$runner\" --trait ${trait} --DGE ${DGE} --cohort ${cohort} --method ${method} --cell_type ${cell_type} --sub_outpath \"${sub_outpath}\" --block ${block} "
          
        sleep 0.2
      fi

    done
  done
done
