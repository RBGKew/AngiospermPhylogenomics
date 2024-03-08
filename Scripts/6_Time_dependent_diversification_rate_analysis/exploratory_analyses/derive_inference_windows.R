library(phytools)

tree <- read.tree("young_tree_smoothing_10_pruned_for_diversification_analyses.tre") # or old_tree_smoothing_10_pruned_for_diversification_analyses.tre

tree_table <- cbind(tree[[1]], tree$edge.length)
tree_table <- tree_table[order(tree_table[,2]),]
tree_table <- cbind(tree_table, c(rep(0, length(tree$tip.label)), branching.times(tree)[-1]))

branch_end_times_ordered <- rev(sort(branching.times(tree)[-1]))
counter <- vector(mode="numeric", length=0)
break_points <- max(node.depth.edgelength(tree))
for (i in 1:length(branch_end_times_ordered)){
counter <- append(counter, 1)
if ((length(counter) > 50) & (branch_end_times_ordered[[i]] < break_points[[length(break_points)]]-5)){
break_points <- append(break_points, branch_end_times_ordered[[i]] - 0.01)
counter <- vector(mode="numeric", length=0)
}
}

