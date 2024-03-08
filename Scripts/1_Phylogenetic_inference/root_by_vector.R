#!/usr/bin/env Rscript
require(ape)

# Positional arguments are: 1- gene trees in a single file; 2- outgroups file, with one sample per line
args <-  commandArgs(trailingOnly=TRUE)
gene_trees_fn <- args[1]
outgroups_fn <- args[2]

root_by_matching_vector <- function(phy, outgroups){
  tips_match <- match(outgroups, phy$tip.label)
  if(all(is.na(tips_match))){
    return(phy)
  }
  first_match <- which(!is.na(tips_match))[1]
  phy_outgroup <- phy$tip.label[tips_match[first_match]]
  root(phy, phy_outgroup, resolve.root = TRUE)
}

gene_trees <- read.tree(gene_trees_fn)
outgroups <- scan(outgroups_fn, character())

rooted_gene_trees <- lapply(gene_trees, FUN = root_by_matching_vector, outgroups = outgroups)

write.tree(rooted_gene_trees, paste0("rooted_",gene_trees_fn))
