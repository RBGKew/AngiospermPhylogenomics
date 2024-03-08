#!/usr/bin/env Rscript

library(ape)
library(adephylo)

identify_outliers <- function(dists){
  dists[dists == 0] <- NA
  dists[dists < 4] <- 1
  overall_mean <- mean(dists, na.rm = TRUE)
  dists_means <- colMeans(dists, na.rm = TRUE)
  pre_outliers <- names(boxplot.stats(dists_means)$out)
  outliers <- pre_outliers[dists_means[pre_outliers] > overall_mean]
  outliers
}

args <- commandArgs(TRUE)
metadata_fn <- args[1]
in_tree_dir <- args[2]
in_tree_suffix <- args[3]
out_tree_dir <- args[4]
log_dir <- args[5]
plot_dir <- args[6]

if(!dir.exists(out_tree_dir)){
  print("Creating output directory")
  dir.create(out_tree_dir)
}

if(is.na(log_dir)){
  log_dir <- in_tree_dir
  print(paste0("Plotting logs in :", in_tree_dir))
}
 
plot_phylo <- !is.na(plot_dir)
if (plot_phylo & !file.exists(plot_dir)){
  plot_phylo <- FALSE
  print("Directory for plots not found. Skipping plots.")
}

metadata <- read.csv(metadata_fn)
metadata$translated_label <- paste(substr(metadata$Order, 1, 3), metadata$Family, metadata$Genus, metadata$Species, metadata$hyb_genes, metadata$hyb_length, metadata$label, sep = "_")
print("Metadata file loaded.")

genetree_fns <- dir(in_tree_dir, pattern = in_tree_suffix, full.names = TRUE)
print(genetree_fns)

for (i in 1:length(genetree_fns)){
  genetree_fn <- genetree_fns[i]
  gene <- gsub(in_tree_suffix, "", basename(genetree_fn))
  print(paste0("Analysing gene ", gene))
  genetree <- read.tree(genetree_fn)
  genetree_metadata <- metadata[metadata$label %in% genetree$tip.label, ]
  orders_count <- table(genetree_metadata$Order)
  orders_to_test <- names(orders_count)[orders_count > 2]
  tree_dists <- as.matrix(distTips(genetree, method = "nNodes"))
  
  outliers <- NULL
  for (ord in orders_to_test){
    ord_tips <- genetree_metadata$label[genetree_metadata$Order == ord]
    if(is.monophyletic(genetree, ord_tips, reroot = TRUE)){
      print(paste0("Skipping monophyletic ", ord))
      next
    }
    ord_dists <- tree_dists[ord_tips, ord_tips]
    outliers <- append(outliers, identify_outliers(ord_dists))
  }
  outliers_fn <- paste0(log_dir, "/", gene, "_backbone_outliers.txt")
  cat(outliers, file = outliers_fn, sep = "\n")
  
  genetree_out_fn <- paste0(out_tree_dir, "/", gene, "_backbone.tre")
  if(length(outliers) > 0){
    genetree_trimmed <- drop.tip(genetree, outliers)
    write.tree(genetree_trimmed, genetree_out_fn)
  } else {
    write.tree(genetree, genetree_out_fn)
  }
  
  if(plot_phylo){
    tip_colour <- rep("black", length(genetree$tip.label))
    tip_colour[match(outliers, genetree$tip.label)] <- "red"
    genetree$tip.label <- paste0(genetree_metadata$translated_label[match(genetree$tip.label, as.character(genetree_metadata$label))])
    
    plot_fn <- paste0(plot_dir, "/", gene, "_backbone_outliers.pdf")
    pdf(plot_fn, 7, 40)
    plot(genetree, cex = 0.2, tip.color = tip_colour)
    nodelabels(genetree$node.label, frame = "none", cex = 0.2, adj = 0.1, col = "red")
    dev.off()
  }
}
