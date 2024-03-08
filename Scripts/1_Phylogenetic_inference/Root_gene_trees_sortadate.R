#!/usr/bin/env Rscript
library(ape)

# Positional arguments are: 1-Species tree filename; 2-gene trees directory; 3(optional)-outgroup of species tree
args <-  commandArgs(trailingOnly=TRUE)
spt_fn <- args[1]
gt_dir <- args[2]

print("Checkpoint 1")

spt <- read.tree(spt_fn)

# If an outgroup is provided, tree is rootted, otherwise just save as ladderized
if (length(args) > 2){
  spt_outgroup <- args[3]
  spt_root <- ladderize(root(spt, spt_outgroup, resolve.root = TRUE), right = FALSE)
} else {
  spt_root <- ladderize(spt, right = FALSE)
}
  write.tree(spt_root, "spt_rooted.tre")

# Create list of tip order. Export is done to be used with phyx's pxrr
is_tip <- spt_root$edge[,2] <= length(spt_root$tip.label)
spt_ordered_tips <- spt_root$tip.label[spt_root$edge[is_tip, 2]]

cat(paste(spt_ordered_tips, sep="", collapse=","), file = "outgroup_string.txt")

# Root trees and save on subdirectory "rooted" in gene trees directory
gt_fns <- dir(gt_dir, include.dirs = FALSE, no.. = TRUE)
dir.create(paste0(gt_dir,"/rooted"), showWarnings = FALSE)
print(gt_fns)

print("Checkpoint 2")


for (gt_fn in gt_fns){
  print("Checkpoint 3")
  gt <- read.tree(paste0(gt_dir,"/",gt_fn))
  
  tips_match <- match(spt_ordered_tips, gt$tip.label)
  first_match <- which(!is.na(tips_match))[1]
  gt_outgroup <- gt$tip.label[tips_match[first_match]]
  gt_rooted <- root(gt, gt_outgroup, resolve.root = TRUE)
  message(paste0("Rooting ", gt_fn, " in ", gt_outgroup))
  write.tree(gt_rooted, paste0(gt_dir,"/rooted/",gt_fn))
}
