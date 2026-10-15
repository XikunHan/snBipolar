# Single-cell multi-cohort dissection of bipolar disorder


## Data
The source data and processed datasets described in this manuscript are available via the [PsychENCODE Knowledge Portal](https://psychencode.synapse.org/), a platform for accessing data, analyses, and tools generated through grants funded by the National Institute of Mental Health’s PsychENCODE program. Access and reuse are subject to the portal's [data-access and attribution requirements](https://psychencode.synapse.org/DataAccess).

The data and code described in this manuscript are available at: [GitHub](https://github.com/XikunHan/snBipolar) and [Synapse](https://www.synapse.org/#!Synapse:syn51534507).


## Repository structure

```text
snBipolar/
├── README.md
├── cell_type_annotation/
│   └── scBD_cell_types.Rmd
├── spatial_registration/
│   ├── sbatch_sp_registration.R
│   └── sbatch_sp_registration.sh
├── differential_expression/
│   ├── create_cell_type_file.R
│   ├── DGE_run_submit.sh
│   ├── DGE_run_each_cell_type.R
│   ├── sc_helper.R
│   ├── DEG_meta_analysis.R
│   └── scBD_DEGs.Rmd
├── eqtl/
│   ├── sceQTL_run_submit.sh
│   └── sceQTL_run_each_cell_type.R
├── mendelian_randomization/
│   ├── clump_gwas.py
│   ├── scMR_run_submit.sh
│   ├── scMR_run_each.sh
│   ├── scMR_run_each.R
│   └── scBD_scMR.Rmd
├── smr/
│   ├── SMR_run_submit.sh
│   └── SMR_run_each.sh
├── cell_cell_communication/
│   ├── cellchat_run_submit.sh
│   └── cellchat_run_each.R
├── gwas_integration/
│   ├── scBD_link_GWAS_loci.Rmd
│   └── f_manhattan.R
└── data/
    ├── DEG_MASTRE_McL_gene_test.csv
    ├── DEG_MASTRE_Sinai_gene_test.csv
    ├── gene_list_50_test.csv
    ├── meta_analysis_IVW_MASTRE.csv
    └── cell_types_25.txt
```


## Running the scripts

Run commands from the `snBipolar/` repository root.

Shell scripts pass the repository path to analysis jobs through `SNBD_ROOT`. R scripts record this path before calling `setwd()`.

Activate an R environment with the packages listed in the script headers before running a launcher.

```bash
# Example commands for individual analysis stages.
Rscript differential_expression/create_cell_type_file.R
bash differential_expression/DGE_run_submit.sh
bash eqtl/sceQTL_run_submit.sh
bash mendelian_randomization/scMR_run_submit.sh
bash smr/SMR_run_submit.sh
bash spatial_registration/sbatch_sp_registration.sh
bash cell_cell_communication/cellchat_run_submit.sh
```

The launchers submit each analysis independently. Prepare the required input files before starting the corresponding stage.

Full expression objects, genotype matrices, GWAS summary statistics, LD references and publication tables are not included in this repository. Input paths are listed in the script headers.

For an R Markdown notebook, start R in the repository root and render after preparing its inputs:

```r
Sys.setenv(SNBD_ROOT = normalizePath("."))
rmarkdown::render(
  "cell_type_annotation/scBD_cell_types.Rmd",
  knit_root_dir = getwd()
)
```

## Analysis scripts

### Cell-type annotation

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [scBD_cell_types.Rmd](cell_type_annotation/scBD_cell_types.Rmd) | Plot an ACTIONet embedding and selected marker genes from `publication/data/sce_snBD_McLean_MtSinai.rds` and `celltype_stats.tsv`. | UMAP and violin figures. |

### Spatial registration

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [sbatch_sp_registration.sh](spatial_registration/sbatch_sp_registration.sh) | Submit the McLean and MtSinai cohort jobs with editable SLURM resources. | Cohort directories and logs beneath `out/sp_registration/`. |
| [sbatch_sp_registration.R](spatial_registration/sbatch_sp_registration.R) | Read `data/sce_snBD_McLean_MtSinai.rds`, fetch the spatialLIBD reference and call `registration_wrapper()` for the requested cohort. | Gene annotation CSV and `sce_modeling_results_<cohort>.rds`. |

### Differential expression

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [create_cell_type_file.R](differential_expression/create_cell_type_file.R) | Split the combined SCE by `Celltype_inferred` and cohort. | `data/cell_type_data_McLean/` and `data/cell_type_data_Sinai/`, each containing `snBD_sce_format_<cell_type>.rds`. |
| [DGE_run_submit.sh](differential_expression/DGE_run_submit.sh) | Submit the selected model/cohort/cell-type/block jobs. | Model-specific result directories under `out/DGE/test/` and job logs. |
| [DGE_run_each_cell_type.R](differential_expression/DGE_run_each_cell_type.R) | Load a split SCE; apply cell-level QC and cohort-specific covariates; estimate three RUVr factors; run edgeR, NEBULA or MAST-RE blocks through `sc_helper.R`. | `<DGE>_DGE__<trait>__<cell_type>[__<block>].csv`. |
| [sc_helper.R](differential_expression/sc_helper.R) | Define `deg.edger()`, `deg.ruvseq()`, `deg.nebula()`, `deg.MAST()`, `se_make_pseudobulk()`, `calculate_qc_metrics()`, `make_pseudobulk()` and `msub()`. | Source this file; it is not a standalone analysis. |
| [DEG_meta_analysis.R](differential_expression/DEG_meta_analysis.R) | Combine cohort effect estimates by inverse-variance weighting and cohort P values using Fisher's method. Uses the two example CSVs by default. | The result is stored in `dd`. Uncomment the final `fwrite()` call to write a CSV. |
| [scBD_DEGs.Rmd](differential_expression/scBD_DEGs.Rmd) | Read the DEG, bulk-expression and enrichment tables from `publication/data/DEG/` and plot overlaps and GO results. | Heatmaps and GO figures. The GO section requires `publication/data/f_reduce_go_term.R`, which is not included. |

For the small meta-analysis example, install the packages listed in its header, then run:

```r
# Current working directory must be the repository root.
source("differential_expression/DEG_meta_analysis.R")
head(dd)
```

For the full analysis, replace the two example input paths with the merged cohort results. The DGE launcher is configured for McLean, Endo and MAST block 1; set the cohort, cell-type and block selections for the analysis being run.


### cis-eQTL mapping

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [sceQTL_run_submit.sh](eqtl/sceQTL_run_submit.sh) | Submit one job for each label in `data/cell_types_25.txt`. | Logs and cell-type results under `out/sceQTL430/result/`. |
| [sceQTL_run_each_cell_type.R](eqtl/sceQTL_run_each_cell_type.R) | Decode genotypes using genio and run linear MatrixEQTL with a 1-Mb cis window. | `genotype_<cell_type>.txt` and `out__<cell_type>__cis_eqtls.txt.gz`. |


### Mendelian randomization

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [clump_gwas.py](mendelian_randomization/clump_gwas.py) | PLINK1.9 LD clumping and selection of summary-statistic rows. Set column names and thresholds through command-line arguments. | PLINK files and `indep_GWhits_*` files. Use `--help` for arguments. `PLINK` can specify the executable. |
| [scMR_run_submit.sh](mendelian_randomization/scMR_run_submit.sh) | Submit one worker per selected cell type. | Logs under `out/sceQTL430/scMR/IV/clump/logs/`. |
| [scMR_run_each.sh](mendelian_randomization/scMR_run_each.sh) | Run BD and SZ MR from clumped instruments. Gene splitting, clumping and merging are optional preparation stages and are disabled by default. | Prepared instruments beneath `scMR/IV/clump/` and results in `scMR/result/`.|
| [scMR_run_each.R](mendelian_randomization/scMR_run_each.R) | Format and harmonize clumped eQTL/GWAS statistics and run the Wald-ratio, IVW, weighted-median/mode and MR-Egger methods. | MR estimates, heterogeneity, intercept-test and single-SNP CSVs. |
| [scBD_scMR.Rmd](mendelian_randomization/scBD_scMR.Rmd) | Read prepared merged MR/SMR tables; perform disease enrichment and plot method overlap by cell type/direction. | Disease-enrichment PDF and comparison plots. |


### SMR

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [SMR_run_submit.sh](smr/SMR_run_submit.sh) | Submit one SMR job per cell type. | Logs under `out/sceQTL430/SMR/logs/`. |
| [SMR_run_each.sh](smr/SMR_run_each.sh) | Run the `--smr-multi` commands on prepared BESD files, LD reference and SMR-format GWAS summaries. | `SMR_out/smr_<trait>_<cell_type>_multi.*`. Input preparation and BESD conversion are disabled by default. |


### Cell-cell communication

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [cellchat_run_submit.sh](cell_cell_communication/cellchat_run_submit.sh) | Submit the McLean/Sinai CellChat jobs. | `out/cellchat/<cohort>/major/` and logs. |
| [cellchat_run_each.R](cell_cell_communication/cellchat_run_each.R) | Read `data/synapse_BD/sce_snBD_McLean_MtSinai.rds`; infer BD and CON networks with `CellChatDB.human`; compare interaction counts and weights and plot signaling pathways. | `cell_chat.RData` and `cell_chat_2.RData` contain lists saved with `saveRDS()`; load them with `readRDS()`. |


### GWAS integration

| Script | Purpose and inputs | Outputs |
|---|---|---|
| [scBD_link_GWAS_loci.Rmd](gwas_integration/scBD_link_GWAS_loci.Rmd) | Combine the annotated Manhattan plot and signed-expression/FDR heatmap from prepared `publication/data/GWAS/` inputs. | GWAS-gene figures.  |
| [f_manhattan.R](gwas_integration/f_manhattan.R) | Define the Manhattan plotting function adapted from Manhattan++ and its helper functions. | Returns `list(final.plot, gt, p)`; plots are not saved automatically.  |
