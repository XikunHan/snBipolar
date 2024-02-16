# # XH, 2023

f_reduce_go_term <- function(data, threshold = 0.7, verbose = TRUE) {
  # data: term_id; P
  cat("similarity terms\n")
  
  simMatrix <- rrvgo::calculateSimMatrix(data$term_id,
                                         orgdb="org.Hs.eg.db",
                                         ont="BP",
                                         method="Rel")
  
  
  set.seed(0)
  scores = sapply(data$P, function(pval) max(-log10(pval)))
  names(scores) = data$term_id
  data$scores = scores
  
  reducedTerms <- rrvgo::reduceSimMatrix(simMatrix,
                                         scores,
                                         threshold=threshold,
                                         orgdb="org.Hs.eg.db")
  
  reducedTerms <- reducedTerms[reducedTerms$parentTerm  %in% reducedTerms$parentTerm[duplicated(reducedTerms$parentTerm)], ]
  
  p_reducedTerms <- rrvgo::treemapPlot(reducedTerms)
  
  if(verbose) {
    print(p_reducedTerms)
  }
  return(reducedTerms)
}

