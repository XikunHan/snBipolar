#!/bin/bash
#SBATCH --job-name module
#SBATCH -c 1                             
#SBATCH -p serial_requeue
#SBATCH -t 10:00   
#SBATCH --mem=10M
#SBATCH -o ./log/run_%x_%N_%j.out
#SBATCH -e ./log/run_%x_%N_%j.err


set -eu 



PATHTOSPACK=./Lab/xikun/software/spack
source $PATHTOSPACK/share/spack/setup-env.sh
spack load anaconda3

conda activate ./Lab/xikun/software/hdWGCNA


threads=1
mem=380000 
walltime=10:00:00
partition=serial_requeue




dir_main=./snBD/out/module/all_gene 

for cell_type in Ex In
do

for cohort in McLean
do


sub_outpath=${dir_main}/${cohort}/${cell_type}
log_path=${sub_outpath}/logs


mkdir -p ${sub_outpath}
mkdir -p ${log_path}


sbatch --job-name "Modu_${cell_type}_${cohort}" --export=ALL \
-e ${log_path}/${cell_type}_${cohort}.err -o ${log_path}/${cell_type}_${cohort}.out \
--cpus-per-task=${threads} --time=${walltime} --mem=${mem} --partition=${partition} \
./sbatch_module.R  --cell_type ${cell_type} --dir_main ${dir_main}  --sub_outpath ${sub_outpath} --cohort ${cohort};


sleep 1
echo -e "\n\n"

done
done

