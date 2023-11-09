


# Hypergeometric test function
hypergeo_test <- function(degs, module_genes, all_genes) {
  overlap <- length(intersect(degs, module_genes))
  m <- length(degs)
  k <- length(module_genes)
  N <- all_genes
  p_value <- phyper(overlap-1, k, N-k, m, lower.tail = FALSE)
  return(list(overlap = overlap,total_DEGs = m, module_Genes = k, All_genes = N, p_value = p_value))
}

# Add significance column
add_significance <- function(p) {
  ifelse(p < 0.001, '***', 
         ifelse(p < 0.01, '**', 
                ifelse(p < 0.05, '*', '')))
}


# Rename contrasts
rename_contrast <- function(x) {
  if (grepl("IVM11vsCtrl", x)) {
    gsub("IVM11vsCtrl", "~ 10^-11 ~M", x)
  } else if (grepl("IVM9vsCtrl", x)) {
    gsub("IVM9vsCtrl", "~ 10^-9 ~M", x)
  } else {
    x
  }
}

