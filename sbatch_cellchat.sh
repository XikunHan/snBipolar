#!/bin/bash
#SBATCH --job-name module
#SBATCH -c 1                             
#SBATCH -p serial_requeue
#SBATCH -t 10:00   
#SBATCH --mem=10M
#SBATCH -o ./log/run_%x_%N_%j.out
#SBATCH -e ./log/run_%x_%N_%j.err



set -eu 


module purge
PATHTOSPACK=./Lab/xikun/software/spack
source $PATHTOSPACK/share/spack/setup-env.sh
spack load anaconda3
source activate ./Lab/xikun/software/renv



threads=1
mem=100000 
walltime=4:00:00
partition=serial_requeue


dir_main=./snBD/out/cellchat


for cell_type in major
do

for cohort in McLean Sinai
do

sub_outpath=${dir_main}/${cohort}/${cell_type}
log_path=${sub_outpath}/logs


mkdir -p ${sub_outpath}
mkdir -p ${log_path}

sbatch --job-name "cellchat${cell_type}_${cohort}" --export=ALL \
-e ${log_path}/${cell_type}_${cohort}.err -o ${log_path}/${cell_type}_${cohort}.out \
--cpus-per-task=${threads} --time=${walltime} --mem=${mem} --partition=${partition} \
./sbatch_cellchat.R  --cell_type ${cell_type}  --sub_outpath ${sub_outpath} --cohort ${cohort};


sleep 1
echo -e "\n\n"

done
done

