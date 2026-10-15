################################################################################
# sc_helper.R — DEG helpers for snBD analyses
#
# Purpose
#   Reusable R helpers for differential expression (DE) analyses on single-nucleus
#   RNA-seq data. Implements:
#     • deg.edger()  : edgeR GLM (QL or LRT) on pseudobulk
#     • deg.nebula() : NEBULA GLMM (cell-level, random effects)
#     • deg.MAST()   : MAST (cell-level, with optional random intercept per donor, MAST-RE)
#     • deg.ruvseq(): residual-based RUVr covariates (added to colData)
#     • se_make_pseudobulk(), make_pseudobulk(): donor aggregation
#     • calculate_qc_metrics(): scanpy-style QC fields
#
# Design choices
#   • Pseudobulk (edgeR): aggregate counts by donor.
#   • Cell-level models: NEBULA uses counts; MAST uses log2-transformed expression.
#     Both support donor-level random effects.
#   • Unwanted variation: RUV (K configurable) appended as RUVW_* covariates.
#   • Output columns are method-specific: MAST returns logFC; edgeR/NEBULA
#     provide log2FC fields. Do not assume identical effect units across methods.
#
# Conventions & expectations
#   • Input is a (Single)SummarizedExperiment with assay "counts" unless overridden.
#   • The pathology column has two levels: <control> and <case>.
#   • Covariates must be present in colData(). The functions do not
#     automatically remove covariates or check the design matrix.
#
# Acknowledgements
#   Helpful discussions with Kellis Lab members; related utilities
#   adapted from lab resources: https://github.com/KellisLab

# Date: 2023-2
# R version: 4.2.0 (2022-04-22) on x86_64-pc-linux-gnu
################################################################################

library(data.table)



# 1. edgeR differential expression -------------------------------------
#' Run edgeR on a SummarizedExperiment
#'
#' @param se SummarizedExperiment to run
#' @param pathology Column in colData(se) that contains case and control
#' @param case Case string within pathology column. Example "AD"
#' @param control Control string within pathology column. Example "control"
#' @param covariates Vector of covariates to add to design matrix
#' @param method Method used for edgeR. Either LRT or QL
#' @param assay Assay from se to utilize
#' @param cpm.cutoff Cutoff of CPM to utilize
#' @param cpm.count Number of observations passing CPM cutoff for filter
#' @export
deg.edger <- function(se, pathology, case, control, covariates=c(),
                      method="LRT", assay=NULL, cpm.cutoff=1, cpm.count=3,
                      filter_only_case_control=FALSE) {
    if (filter_only_case_control) {
        se = se[,SummarizedExperiment::colData(se)[[pathology]] %in% c(case, control)]
    }
    
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)$counts
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    if (ncol(X) > 1000) {
        warning("More than 1000 samples are being called to a function under the assumption data passed in is pseudobulk")
    }
    
    cd = SummarizedExperiment::colData(se)
    if (!is.character(cd[[pathology]]) | !is.factor(cd[[pathology]])) {
        warning(paste0("Converting ", pathology, " to factor"))
    }
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    dgel = edgeR::DGEList(X, group=cd[[pathology]], remove.zeros=TRUE)
    to_keep = rowSums(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.count
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    dgel = edgeR::calcNormFactors(dgel, method="TMM")
    design = model.matrix(as.formula(paste0("~0 + ", paste0(c(pathology, covariates), collapse=" + "))),
                          data=cd)
    colnames(design) = make.names(colnames(design)) ### need spaces to be removed
    contrasts = limma::makeContrasts(contrasts=paste0(make.names(paste0(pathology, case)),
                                                      "-",
                                                      make.names(paste0(pathology, control))),
                                     levels=design)
    if (!(method %in% c("QL", "LRT"))) {
        method = "QL"
        warning(paste0("Setting method to ", method))
    }
    if (method == "QL") {
        dgel = edgeR::estimateDisp(dgel, design, robust=TRUE)
        fit <- edgeR::glmQLFit(dgel, design, robust=TRUE)
        res <- edgeR::glmQLFTest(fit, contrast=contrasts)
    } else if (method == "LRT") {
        dgel = edgeR::estimateGLMCommonDisp(dgel, design)
        dgel = edgeR::estimateGLMTagwiseDisp(dgel, design)
        fit = edgeR::glmFit(dgel, design, robust=TRUE)
        res = edgeR::glmLRT(fit, contrast=contrasts)
    }
    tbl = edgeR::topTags(res, n=Inf, sort.by="none")$table
    colnames(tbl) = gsub("^logFC$", "log2FC", colnames(tbl))
    tbl$gene = rownames(tbl)
    tbl$method = method
    tbl$case = case
    tbl$control = control
    return(tbl[order(tbl$FDR),])
}

# 2. Residual-based unwanted-variation factors -------------------------
#' Add RUVseq columns to SummarizedExperiment
#' Based on Boix's RUVSeq workflow
#'
#' @param sce SummarizedExperiment object
#' @param sample Sample on which to pseudobulk
#' @param pathology Main variable of interest
#' @param NRUV Number of residual-based RUVr components requested.
#' @param assay Assay from SummarizedExperiment to use
#' @param cpm.cutoff Cutoff of CPM to utilize
#' @param cpm.count Number of observations passing CPM cutoff for filter
#' @param norm edgeR norm method for calcNormFactors
#' @export
deg.ruvseq <- function(sce, sample, pathology, NRUV=10, assay=NULL, cpm.cutoff=1, cpm.count=3, norm="TMM") {
    pb = se_make_pseudobulk(sce, sample)
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(pb)$counts
    } else {
        X = SummarizedExperiment::assays(pb)[[assay]]
    }
    cd = SummarizedExperiment::colData(pb)
    dgel = edgeR::DGEList(X, remove.zeros=TRUE)
    # print(dgel)
    to_keep = Matrix::rowSums(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.count
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    design = model.matrix(as.formula(paste0("~", pathology)), cbind(dgel$samples, cd))
    # print(design)
    ### Use LRT workflow
    dgel = edgeR::calcNormFactors(dgel, method=norm)
    dgel = edgeR::estimateGLMCommonDisp(dgel, design)
    dgel = edgeR::estimateGLMTagwiseDisp(dgel, design)
    fit1 = edgeR::glmFit(dgel, design)
    res1 = residuals(fit1, type="deviance")
    ruv = RUVSeq::RUVr(dgel$counts,
                       cIdx=as.character(rownames(dgel$counts)),
                       k=NRUV,
                       residuals=res1)
    
    ruvw = ruv$W
    colnames(ruvw)  = paste0("RUV", colnames(ruvw))
    rownames(ruvw) = colnames(ruv$normalizedCounts)
    S = as.character(SummarizedExperiment::colData(sce)[[sample]])
    for (cn in colnames(ruvw)) {
        SummarizedExperiment::colData(sce)[[cn]] = ruvw[S, cn]
    }
    return(sce)
}




# 3. NEBULA mixed-model differential expression ------------------------
#' Run Nebula on a SummarizedExperiment
#'
#' @param sce SummarizedExperiment object
#' @param sample Sample which to model dispersions
#' @param pathology Pathology column to study
#' @param case Case to look at within pathology column
#' @param control Control value to look at within pathology column
#' @param covariates Vector of covariate columns within colData(sce)
#' @param offset Offset from which to weight cells in model
#' @param assay Assay within sce. Default is "counts"
#' @param model Model to pass to nebula::nebula
#' @param filter_only_case_control logical telling whether to filter sce to only case and control values within pathology
#' @param NRUV Number of RUVSeq::RUVr components. 0 means no RUVr components
#' @param cpm.cutoff Cutoff of CPM to utilize for RUV
#' @param cpm.count Number of observations passing CPM cutoff for filter for RUV
#' @export
#' 

deg.nebula <- function(sce, sample, pathology, case, control, covariates=c(),
                       offset="total_counts", assay=NULL, model="NBGMM",
                       filter_only_case_control=TRUE, factorize_pathology=TRUE,
                       cpm.cutoff=1, cpm.count=3, NRUV=10) {
    if (filter_only_case_control) {
        sce = sce[,SummarizedExperiment::colData(sce)[[pathology]] %in% c(case, control)]
    }
    if (NRUV > 0) {
        sce = deg.ruvseq(sce,
                         sample=sample,
                         pathology=pathology,
                         NRUV=NRUV,
                         assay=assay,
                         cpm.cutoff=cpm.cutoff,
                         cpm.count=cpm.count)
        cd = SummarizedExperiment::colData(sce)
        covariates = c(covariates, colnames(cd)[grep("^RUV", colnames(cd))])
    }
    sce[[sample]] = as.factor(as.character(sce[[sample]]))
    sce = sce[,order(sce[[sample]])]
    cd = SummarizedExperiment::colData(sce)
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(sce)$counts
    } else {
        X = SummarizedExperiment::assays(sce)[[assay]]
    }
    if (factorize_pathology) {
        cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    }
    design = model.matrix(as.formula(paste0("~",
                                            paste0(c(pathology, covariates), collapse="+"))),
                          data=cd)
    if (is.null(offset) || !(offset %in% colnames(cd))) {
        warning("Offset not present. Using number of counts")
        offset = Matrix::colSums(X)
    } else {
        print(paste0("Using offset ", offset))
        offset = cd[[offset]]
    }
    neb = nebula::nebula(X,
                         cd[[sample]],
                         design,
                         offset=offset,
                         model=model)

    
    ### Now to parse output
    name.case = paste0(pathology, case)
    name.control = paste0(pathology, control)
    if (paste0("p_", name.case) %in% colnames(neb$summary)) {
        neb$summary$FDR = p.adjust(neb$summary[[paste0("p_", name.case)]], "fdr")
        neb$summary$log2FC = neb$summary[[paste0("logFC_", name.case)]] / log(2)
    } else if (paste0("p_", name.control) %in% colnames(neb$summary)) {
        neb$summary$FDR = p.adjust(neb$summary[[paste0("p_", name.control)]], "fdr")
        neb$summary$log2FC = -1 * neb$summary[[paste0("logFC_", name.control)]] / log(2)
    }
    neb_df = neb$summary
    ovr = neb$overdispersion
    colnames(ovr) = paste0("overdispersion_", colnames(ovr))
    neb_df = cbind(neb_df, ovr)
    neb_df$convergence = neb$convergence
    neb_df$algorithm = neb$algorithm
    neb_df$case = case
    neb_df$control = control
    return(neb_df[order(neb_df$FDR),])
}





###############################################################################
# snBD — Function: deg.MAST()
#
# Purpose:
#   Fit the MAST hurdle model, with optional donor random intercept.
#
# Inputs:
#   se: (Single)SummarizedExperiment; pathology/case/control: contrast labels.
#   sample: donor field; covariates: colData fields; assay: input assay name.
#   detection_rate: include standardized detected-gene count.
#   method="glmer": include a donor random intercept.
#   gene_list: optional gene subset; min_cell: detected-cell threshold (0 disables).
#   ebayes and filter_only_case_control: model and sample-filtering options.
#
# Outputs:
#   Table with gene, P (hurdle test), logFC, ci.hi, ci.lo, FDR, case, control.
#
# Method overview:
#   Calculate normalized expression; construct log2(X+1); build the model;
#   optionally select a gene block; run MAST::zlm and summarize the contrast.
#
###############################################################################


deg.MAST <- function(se,
  pathology,
  sample="SubID",
  case,
  control,
  covariates=c(),
  assay=NULL,
  detection_rate = TRUE,
  method = "glmer", 
  ebayes = FALSE, 
  filter_only_case_control=FALSE,
  gene_list = NULL,
  min_cell = 0) {
  
  if (filter_only_case_control) {
    se = se[,SummarizedExperiment::colData(se)[[pathology]] %in% c(case, control)]
  }
  
# Optional gene filter: expression must be detected in more than min_cell cells.
  if(min_cell > 0) {
    
    d <- scater::nexprs(se, byrow=TRUE)
    se <- se[d > min_cell, ]
    
  }
  
  if (is.null(assay)) {
    X = SummarizedExperiment::assays(se)$counts
  } else {
    X = SummarizedExperiment::assays(se)[[assay]]
  }
  
  # X <- X[which(apply(X, 1, mean) > 0), ]
  
  # TPM 
  library(MAST)
  # devtools::install_github("RGLab/MAST")
  
  X <- scuttle::calculateTPM(X)
  sca <- MAST::FromMatrix(exprsArray = as.matrix(log2(X + 1)), 
    cData = data.frame(wellKey=colnames(se), colData(se)))
  
  colData(sca)$cdr <- c(scale(colSums(assay(sca)>0)))
  
  colData(sca)[[sample]]  <- factor(colData(sca)[[sample]])
  
  if(is.null(covariates)) {
    formula <- paste0(" ~ ", pathology)
  } else {
    formula = paste0(" ~  ", paste0(c(pathology, covariates), collapse=" + "))
  }
  
  
  # covariates
  if(detection_rate) {
    formula <- paste0(formula, " + cdr") # Include the detection-rate covariate when requested.
  }
  
  # sca$outcome <- factor(sca$outcome)
  # sca$outcome <- relevel(sca$outcome, ref = control)
  # 
  print(sca)
  # print(sca@colData)
  print(str(sca@colData))
  
  if(is.null(gene_list)) {
    
    if(method == "glmer") {
      # Runtime increases with the number of cells. 
      print(sprintf("Running %s", method))
      formula <- paste0(formula, sprintf("+ (1|%s)", sample))
      print(as.formula(formula)) 
      
      library(lme4)
      zlmdata <- MAST::zlm(as.formula(formula), sca, method = method, ebayes = ebayes, strictConvergence = FALSE)
      
      
    } else if(method == "general") {
      
      print(sprintf("Running %s", method))
      print(as.formula(formula)) 
      
      zlmdata <- MAST::zlm(as.formula(formula), sca)
    }
    
  } else {
    # Fit the random-effects model to the selected gene block.
    if(method == "glmer") {
      print(sprintf("Running %s", method))
      formula <- paste0(formula, sprintf("+ (1|%s)", sample))
      print(as.formula(formula)) 
      
      sca <- sca[rownames(sca) %in% gene_list, ]
      
      print(sca)
      
      library(lme4)
      zlmdata <- MAST::zlm(as.formula(formula), sca, method = method, ebayes = ebayes, strictConvergence = FALSE)
      
    }
    
  }

  
  v_test_var <- paste0(pathology, case)
  
  summary_res <- summary(zlmdata, doLRT= v_test_var) 
  
  summaryDt <- summary_res$datatable
  
  fcHurdle <- merge(summaryDt[contrast== v_test_var  & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
    summaryDt[contrast== v_test_var & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  names(fcHurdle) <- c("gene", "P", "logFC", "ci.hi", "ci.lo", "FDR")
  fcHurdle <-  fcHurdle[order(fcHurdle$FDR),]
  fcHurdle$case = case
  fcHurdle$control = control
  return(fcHurdle)
}








# 5. SummarizedExperiment pseudobulk aggregation -----------------------
#' Pseudobulk a SummarizedExperiment
#'
#' Determines constant columns across dataframe
#' and summarizes counts
#'
#' @param se SummarizedExperiment
#' @param on Column in colData() to split upon
#' @export
se_make_pseudobulk <- function(se, on) {
    cd = SummarizedExperiment::colData(se)
    rd = SummarizedExperiment::rowData(se)
    A = SummarizedExperiment::assays(se)
    to_P = Matrix::Matrix(make_pseudobulk(cd[[on]]))
    P = lapply(A, function(mat) {
        mat %*% to_P
    })
    cd[[on]] = factor(as.character(cd[[on]]), levels=colnames(to_P))
    pcd = data.frame(row.names=levels(cd[[on]]))
    for (cn in colnames(cd)) {
        ### for each "on", ensure all equal
        u = lapply(split(cd[[cn]], cd[[on]]), unique)
        good = all(sapply(u, function(x) { length(x) == 1 }))
        if (good) {
            pcd[[cn]] = sapply(u, identity)[rownames(pcd)]
            if (is.factor(cd[[cn]])) {
                pcd[[cn]] = as.factor(as.character(pcd[[cn]]))
            }
        }
    }
    return(SummarizedExperiment::SummarizedExperiment(assays=P,
                                                      rowData=rd,
                                                      colData=pcd))
}

# 6. Expression QC metrics ---------------------------------------------
#' Calculate qc metrics similar to Scanpy
#'
#' @param se SummarizedExperiment
#' @param assay Assay to calculate upon
#' @param log1p logical whether to use log statistics as well
#' @export
calculate_qc_metrics <- function(se, assay=NULL, log1p=TRUE) {
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)[[1]]
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    rs = Matrix::rowSums(X)
    SummarizedExperiment::rowData(se)$total_counts = rs
    SummarizedExperiment::rowData(se)$n_cells_by_counts = Matrix::rowSums(X > 0)
    cs = Matrix::colSums(X)
    csz = Matrix::colSums(X > 0)
    SummarizedExperiment::colData(se)$total_counts = cs
    SummarizedExperiment::colData(se)$n_genes_by_counts = csz
    SummarizedExperiment::colData(se)$n_genes = csz
    if (log1p) {
        SummarizedExperiment::rowData(se)$log1p_total_counts = base::log1p(rs)
        SummarizedExperiment::colData(se)$log1p_total_counts = base::log1p(cs)
        SummarizedExperiment::colData(se)$log1p_n_genes_by_counts = base::log1p(csz)
    }
    return(se)
}






# 7. Cell-to-donor membership matrix -----------------------------------
#' Collapse a matrix into pseudobulk
#'
#' @param x Vector of cell-level group labels (not a count matrix)
#' @param u unique levels to set
#' @return Cell-by-group membership matrix used for subsequent multiplication
#' @export
make_pseudobulk <- function(x, u=NULL) {
    if (is.factor(x)) {
        u = levels(x)
    } else if (is.null(u)) {
        u = unique(x)
    }
    u = u[!is.na(u)]
    names(u) = u
    return(sapply(u, function(y) {
        1 * (x == y)
    }))
}




###############################################################################
# snBD — Function: msub()
#
# Purpose:
#   Apply successive gsub replacements to a character vector.
#
# Inputs:
#   pattern/replacement: ordered replacements; x: character vector.
#   Options may be scalars or vectors.
#
# Outputs:
#   Character vector after each replacement has been applied in order.
#
###############################################################################


msub <- function(pattern, replacement, x, ignore.case=FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
    # Helper: select element i for a vector option, or reuse its scalar value.
    lgt1i <- function(item, i) {
        if (length(item) > 1) {
            item[i]
        } else {
            item
        }
    }
    for (i in seq_along(pattern)) {
        x = gsub(pattern[i],
                 replacement[i],
                 x,
                 ignore.case=lgt1i(ignore.case, i),
                 perl=lgt1i(perl, i),
                 fixed=lgt1i(fixed, i),
                 useBytes=lgt1i(useBytes, i))
    }
    return(x)
}



