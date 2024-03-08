library(ape)

# Positional arguments are: [1] gene alignment directory [2] gene alignment pattern [3] gene trees directory [4] gene trees pattern [5] destination
args <-  commandArgs(trailingOnly=TRUE)
aligns_dir <- args[1]
aligns_pattern <- args[2]
trees_dir <- args[3]
trees_pattern <- args[4]
trees_destination <- args[5]
trees_out_suffix <- args[6]

aligns_fns <- dir(aligns_dir, pattern = aligns_pattern, include.dirs = FALSE, no.. = TRUE)

for (aln_fn in aligns_fns){
  gene <- gsub(aligns_pattern, "", aln_fn)
  tree_output <- paste0(trees_destination, "/", gene, trees_out_suffix)
  if (file.exists(tree_output)){
    next
  }
  aln <- read.FASTA(paste0(aligns_dir, "/", aln_fn), type = "DNA")
  tree_fn <- paste0(trees_dir, "/", gene, trees_pattern)
  tree <- read.tree(tree_fn)
  tips_excess <- tree$tip.label[!tree$tip.label %in% names(aln)]
  if (length(tips_excess) > 0) {
    tree <- drop.tip(tree, tip = tips_excess)
  }
  write.tree(tree, tree_output)
}
