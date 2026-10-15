#!/bin/bash

###############################################################################
# snBD — Prepare clumped eQTL instruments and run scMR
#
# Purpose:
#   Prepare clumped instruments for one cell type and run MR for BD and SZ.
#
# Inputs:
#   Environment variables cell_type and dir_main set by the launcher.
#   eQTL results, PLINK1.9, European LD reference and BD/SZ GWAS files.
#
# Outputs:
#   out/sceQTL430/scMR/IV/clump/ and scMR/result/ files.
#
# Method overview:
#   Optional preparation: select eQTL FDR<0.05, split by gene, clump at
#   r2=0.01 within 1000 kb with p1=p2=1, merge instruments. Then run BD and SZ MR.
#
# Dependencies:
#   Bash, R/Rscript, Python, PLINK1.9, data.table, TwoSampleMR.
#
# Usage:
#   Called by mendelian_randomization/scMR_run_submit.sh.
#   Activate an R environment with the required packages before submitting.
#
###############################################################################

# Pass the repository root to the analysis jobs.
# 1. Repository root and local environment -----------------------------
export SNBD_ROOT="${SNBD_ROOT:-$PWD}"

set -eu 

# module load R/4.2.0


echo ${cell_type}
export PLINK="${PLINK:-${SNBD_ROOT}/software/plink_linux_x86_64_20201019_v1.9/plink}"
cd "${dir_main}"


mkdir -p  scMR/IV/clump/${cell_type}/data
mkdir -p  scMR/IV/clump/${cell_type}/result
mkdir -p  scMR/IV/clump/merge


if [[ 0 -gt 1 ]]
then


# run the clump methods
if [[ 2 -gt 1 ]]
then

# read data in R, to clean the data
R --vanilla <<EOF
library(data.table)
cell_type <- "${cell_type}"
print(cell_type)

df <- fread(paste0("result/out__", cell_type,"__cis_eqtls.txt.gz"))
df

df[, table(FDR < 0.05)]

df <- df[FDR < 0.05]
df

print(length(unique(df\$gene)))

for(i in seq_along(unique(df\$gene))) {
  
  gene_i <- unique(df\$gene)[i]
  print(paste0(i, "_", gene_i))
  
  df_s <- df[gene == gene_i]
  df_s
  
  file <- paste0("scMR/IV/clump/", cell_type, "/data/", cell_type, "__", gene_i, ".txt")
  fwrite(df_s, file, sep = "\t")
}


EOF



# module load GCCcore/8.2.0
# module load Python/3.7.2

# export PATH="./software/anaconda3/bin:$PATH"

echo "Run clumping via PLINK"

files=`ls scMR/IV/clump/${cell_type}/data`

for file in ${files}
do 

gwas_file=scMR/IV/clump/${cell_type}/data/${file}

python "${SNBD_ROOT}/mendelian_randomization/clump_gwas.py" \
--dir_path scMR/IV/clump/${cell_type}/result \
--LDref_bfile ${SNBD_ROOT}/data/1kg_LD_ref/from_cause/EUR \
--gwas_summary_file ${gwas_file} \
--trait ${file} \
--clump_snp snps \
--p1thres 1 \
--p2thres 1 \
--r2level 0.01 \
--clump_field pvalue

done




fi


# merge the files 

if [[ 2 -gt 1 ]]
then

R --vanilla <<EOF

library(data.table)
cell_type <- "${cell_type}"
print(cell_type)

v_clump <- list.files(sprintf("scMR/IV/clump/%s/result", cell_type), pattern = "_1000kb.r2_0.01_p1_1_p2_1")
print(paste0("clump data: ", length(v_clump)))

# Function: f_load_files()
# Purpose: combine clumped gene tables and record each source filename.
# Inputs: path (folder), pattern (filename regular expression).
# Returns: combined data.table; stops when no matching files exist.
f_load_files <- function(path = ".", pattern = NULL) {
  v_files <- list.files(path = path, pattern = pattern)
  print(paste0("file N = ", length(v_files)))
  df_res <- NULL
  if (!length(v_files)) stop("No clumped files found; run instrument preparation first.")
  pb <- txtProgressBar(min = 0, max =length(v_files), style = 3, file = stderr())
  for(j in seq_along(v_files)) {
    setTxtProgressBar(pb = pb, value = j)
    df_one <- data.table::fread(file.path(path, v_files[j]), fill=TRUE)
    df_res <- rbind(df_res, cbind(df_one, file = v_files[j]),  fill=TRUE)
  }
  return(df_res)
}


df_clump <- f_load_files(path = sprintf("scMR/IV/clump/%s/result", cell_type), pattern = "_1000kb.r2_0.01_p1_1_p2_1")
print(paste0("clump results: ", length(unique(df_clump\$file))))
df_clump[, cell_type := cell_type]
fwrite(df_clump, sprintf("scMR/IV/clump/merge/df_merge_clump_%s.txt", cell_type), sep = "\t")

EOF


fi


fi




# module load R/4.0.2

if [[ 2 -gt 1 ]]
then


for trait in BD   SZ
do 
Rscript --slave "${SNBD_ROOT}/mendelian_randomization/scMR_run_each.R" --trait ${trait} --cell_type ${cell_type} --sub_outpath "${dir_main}/scMR/result"

done 

fi
