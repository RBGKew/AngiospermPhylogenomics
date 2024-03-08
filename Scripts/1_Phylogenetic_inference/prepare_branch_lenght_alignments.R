#! /mnt/shared/scratch/azuntini/apps/conda/bin/Rscript

options(show.error.locations = TRUE)
require(optparse)

option_list = list(
  make_option(c("-f", "--gg_file"), action="store", default=NA, type='character',
              help="Path to gg_file"),
  make_option(c("-g", "--gene_thresholds"), action="store", default=NA, type='character',
              help="Gene thresholds. If multiple, provide as '25,50,75,..."),
  make_option(c("-s", "--sp_tree"), action="store", default=NA, type='character',
              help="just a variable named b"),
  make_option(c("-t", "--gene_trees_pattern"), action="store", default=NA, type='character',
              help="Gene tree ending pattern"),
  make_option(c("-a", "--gene_aligns_path"), action="store", default=NA, type='character',
              help="Path to gene alignments"),
  make_option(c("-b", "--gene_aligns_pattern"), action="store", default=NA, type='character',
              help="Gene alignments ending pattern"),
  make_option(c("-r", "--relative_bp"), action="store_true", default=FALSE,
              help="Sort gg_file by relative bipartitions [default %default]"),
  make_option(c("-c", "--constraint_to_smaller"), action="store_true", default=FALSE,
              help="Constraint taxon set to the smallest threshold [default %default]"),
  make_option(c("-l", "--remove_node_labels"), action="store_false", default=TRUE,
              help="Remove node labels [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option = FALSE))

prepare_branch_length_estimation_from_align <- function(gg_file,
                                             gene_thresholds = c(25, 50, 75, 100, 200, 353),
                                             species_tree_path,
                                             gene_tree_pattern,
                                             gene_aligns_path,
                                             gene_aligns_pattern,
                                             sort_by_relative_bp = FALSE,
                                             constraint_to_smaller = FALSE,
                                             remove_node_labels = TRUE){

  print(gg_file)
  library(ape)

  gg <- read.delim(gg_file, sep=" ", header = TRUE)
  if(sort_by_relative_bp){
    gg <- gg[order(gg$bipartition_rel, decreasing = TRUE), ]
  }

  max_genes <- dim(gg)[1]
  gene_thresholds[gene_thresholds > max_genes] <- max_genes
  gene_thresholds <- sort(unique(gene_thresholds))
  max_threshold <- max(gene_thresholds)
  
  species_tree <- read.tree(species_tree_path)
  if(remove_node_labels){
    species_tree$node.label <- NULL
  }
  
  genes <- gsub(gene_tree_pattern, "", gg$name[1:max_threshold])
  aligns_fn <- paste0(gene_aligns_path, "/", genes, gene_aligns_pattern)
  aligns_list <- lapply(aligns_fn, read.FASTA)
  message(paste0(max_threshold, " alignments loaded."))
  
  aligns_th_tips_list <- list()
  
  for (i in 1:length(gene_thresholds)){
    th <- gene_thresholds[i]
    aligns_list_th <- aligns_list[1:th]
    aligns_th_tips_list[[i]] <- unique(unlist(lapply(aligns_list_th, names)))
    
    sub_name <- paste0("top_", th)
    if(sort_by_relative_bp){
      sub_name <- paste0(sub_name, "_rel")
    }
    
    if(constraint_to_smaller){
      align_th_tips <- aligns_th_tips_list[[1]]
      sub_name <- paste0(sub_name, "_small")
    } else {
      align_th_tips <- aligns_th_tips_list[[i]]
    }
    
    th_missing_tips <- species_tree$tip.label[!species_tree$tip.label %in% align_th_tips]
    species_tree_th <- drop.tip(species_tree, th_missing_tips)
    species_tree_th$edge.length <- rep(0.1, length(species_tree_th$edge.length))
    
    if(dir.exists(sub_name)){
      unlink(sub_name, TRUE)
    }
    dir.create(paste0(sub_name, "/aligns"), recursive = TRUE)
    write.tree(species_tree_th, paste0(sub_name,"/", sub_name, "_", "constraint.nwk"))
    
    if(constraint_to_smaller){
      for (j in 1:th){
        align_in <- aligns_list_th[[j]]
        align_out <- align_in[names(align_in) %in% align_th_tips]
        write.FASTA(align_out, paste0(sub_name, "/aligns/", genes[j], ".fasta"))
      }
    } else {
      for (j in 1:th){
        write.FASTA(aligns_list_th[[j]], paste0(sub_name, "/aligns/", genes[j], ".fasta"))
      }
    }
  message(paste0("Threshold ", th, " prepared."))
  }
}  

prepare_branch_length_estimation_from_align(
  opt$gg_file,
  as.numeric(unlist(strsplit(opt$gene_thresholds, ","))),
  opt$sp_tree,
  opt$gene_trees_pattern,
  opt$gene_aligns_path,
  opt$gene_aligns_pattern,
  opt$relative_bp,
  opt$constraint_to_smaller,
  opt$remove_node_labels
  )
