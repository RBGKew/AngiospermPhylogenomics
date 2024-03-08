library(phytools)
library(BAMMtools)

tree <- read.tree("young_tree_smoothing_10_pruned_for_diversification_analyses.tre")
setBAMMpriors(tree, total.taxa = 360000, traits = NULL, outfile = "myPriors.txt", Nmax = 1000, suppressWarning = FALSE)
