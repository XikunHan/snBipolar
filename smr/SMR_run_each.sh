#!/bin/bash

###############################################################################
# snBD — SMR analysis of cell-type eQTLs and BD/SZ GWAS
#
# Purpose:
#   Run the --smr-multi analysis from existing BESD files.
#
# Inputs:
#   Environment variables cell_type, sub_outpath and threads from the launcher.
#   Existing SMR_format/mybesd_<cell_type>.besd/.esi/.epi; LD reference;
#   BD/SZ GWAS files in SMR-compatible summary-statistic format.
#
# Outputs:
#   SMR_out/smr_<trait>_<cell_type>_multi.* under the worker output directory.
#
# Method overview:
#   Input preparation and BESD conversion are disabled. Run --smr-multi
#   separately for BD and SZ using the prepared BESD files.
#
# Dependencies:
#   Bash, SMR executable, PLINK-format LD reference; optional preparation uses R/data.table.
#
# Usage:
#   From repository root: bash smr/SMR_run_submit.sh
#
###############################################################################

# Pass the repository root to the analysis jobs.
# 1. Repository root and local environment -----------------------------
export SNBD_ROOT="${SNBD_ROOT:-$PWD}"

set -eu 

# module load R/4.2.0


echo ${cell_type}

mkdir -p  "${sub_outpath}/data"
mkdir -p  "${sub_outpath}/SMR_format"
mkdir -p  "${sub_outpath}/SMR_out"


cd "${sub_outpath}"

SMR="${SNBD_ROOT}/software/smr-1.3.1-linux-x86_64/smr-1.3.1"
bfile="${SNBD_ROOT}/data/1kg_LD_ref/g1000_eur"


if [[ 0 -gt 1  ]]
then

# read data in R, to clean the data
R --vanilla <<EOF
# cell_type <- "Exc-L6b"

library(data.table)
cell_type <- "${cell_type}"
print(cell_type)


df_freq <- fread("${SNBD_ROOT}/data/1kg_LD_ref/from_cause/g1000_eur_freq.frq")
df_freq

df <- fread(paste0("${SNBD_ROOT}/out/sceQTL430/result/out__",cell_type,"__cis_eqtls.txt.gz"))
df

df[, id := paste0(gene, "_", snps)]

df <- df[order(pvalue)]
df
print(length(df[duplicated(id)]))

df <- df[!duplicated(id)]

# for the eQTL data
df_s <- df[, list(SNP = snps, gene, beta, \`t-stat\` = statistic,  \`p-value\` = pvalue, FDR)]

summary(df_s)
df_s <- df_s[complete.cases(df_s), ]

file <- sprintf("data/df_matrixeqtl__%s.txt", cell_type)
fwrite(df_s, file, sep = "\t")

df_esi <- df[, list(CHR, SNP = snps, tmp = 0, BP, A1 = REF, A2 = ALT)]
df_esi

df_m <- merge(df_esi, df_freq, by = "SNP")

df_m[, freq := ifelse(A1.x == A1.y & A2.x == A2.y, MAF, 1-MAF)]
df_m
df_s <- df_m[, list(CHR = CHR.x, SNP, tmp, BP, A1 = A1.x, A2 = A2.x, freq)]
df_s

df_s <- df_s[!duplicated(SNP)]

summary(df_s)
df_s <- df_s[complete.cases(df_s), ]
fwrite(df_s, sprintf("data/df_matrixeqtl__%s.esi", cell_type), sep = "\t")

df[, pos := as.integer(median(BP)), by = "gene"]
df

df_s <- df[, list(CHR = CHR, 
                  prob = gene,
                  temp = 0,
                  pos = BP,
                  gene = gene,  
                  dir = "+"
)]

df_s <- df_s[!duplicated(gene)]
df_s
summary(df_s)
df_s <- df_s[complete.cases(df_s), ]

fwrite(df_s, sprintf("data/df_matrixeqtl__%s.epi", cell_type), sep = "\t")

EOF
fi


# BESD conversion requires MatrixEQTL-format statistics, allele frequencies
# and gene/probe annotations. Check the annotations before enabling this stage.
if [[ 0 -gt 1  ]]
then

mateQTL="${sub_outpath}/data/df_matrixeqtl__${cell_type}.txt"

"${SMR}" --eqtl-summary "${mateQTL}" --matrix-eqtl-format --make-besd --out "${sub_outpath}/SMR_format/mybesd_${cell_type}"


# update 
"${SMR}" --beqtl-summary "${sub_outpath}/SMR_format/mybesd_${cell_type}" --update-esi "${sub_outpath}/data/df_matrixeqtl__${cell_type}.esi"
"${SMR}" --beqtl-summary "${sub_outpath}/SMR_format/mybesd_${cell_type}" --update-epi "${sub_outpath}/data/df_matrixeqtl__${cell_type}.epi"


fi




if [[ 2 -gt 1  ]]
then


# trait=BD
for trait in BD SZ
do

if [[ $trait == "BD" ]]
then
gwas_summary="${SNBD_ROOT}/out/scMR/gwas/pgc_bip2021_all.txt"

elif [[ $trait == "SZ" ]]
then

gwas_summary="${SNBD_ROOT}/out/scMR/gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.txt"
fi 


"${SMR}" --bfile "${bfile}" \
--gwas-summary "${gwas_summary}" \
--smr-multi \
--beqtl-summary "${sub_outpath}/SMR_format/mybesd_${cell_type}" \
--out "${sub_outpath}/SMR_out/smr_${trait}_${cell_type}_multi" \
--thread-num ${threads}

done


fi
