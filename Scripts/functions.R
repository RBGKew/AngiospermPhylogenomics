#!/usr/bin/env Rscript
#
# functions.R 
#
# functions.R congregates functions to manipulate metadata, datasets, alignments and
# trees used in the analyses of the Big Tree
# 
# 
# #### Version
# 
# - version = 1.0
# - date: 2022.01.18
#
#
# #### Author
# 
# - Alexandre R. Zuntini
# - a.zuntini@kew.org
# 

#### Save as object or file ####

save_object <- function(obj, out_file, out_sep = ","){
  if(is.null(out_file)){
    return(obj)
  } else {
    obj_class <- class(obj)
    if(obj_class == "data.frame"){
      write.table(obj, file = out_file, sep = out_sep, row.names = FALSE, quote = FALSE)
    } else if(obj_class == "character"){
      cat(obj, file = out_file, sep = "\n")
    } else if(obj_class %in% c("Phylo", "phylo", "multiPhylo")){
      write.tree(obj, file = out_file)
    } else {
      print("Unknown format")
    }
  }
}

add_to_factor <- function(fact, value, index){
  if(!value %in% levels(fact)){
    levels(fact) <- c(levels(fact), value)
  }
  fact[index] <- value
  return(fact)
}

is_blank <- function(x){
  return(is.na(x) | x == "")
}


#### Export PAFTOL metadata ####

# 
export_path <- "DB_exports/"
filtering_field <- "idSequencing"

export_metadata <- function(export_path = "DB_exports/",
                            filter_ids = NULL,
                            filtering_field = "idSequencing",
                            fields = NULL,
                            fix_incertae_sedis = TRUE,
                            add_clades = TRUE,
                            recovery_file = NULL,
                            correct_orders = "Supporting_files/new_orders.csv",
                            correct_metadata = NULL,
                            out_file = NULL){
  
  library(stringr)
  
  export_files <- dir(export_path, "_paftol_export.csv", full.names = TRUE)
  latest_export <- export_files[length(export_files)]
  DF <- read.csv(latest_export)
  
  # Clean dataset
  DF <- DF[!is.na(DF$idPaftol), ]
  DF <- DF[DF$Library_Status != "N.A.", ]
  DF <- DF[!DF$Action %in% c("Excluded", "On standby", "Replaced"), ]
  
  # Fix Incertae Sedis
  if(fix_incertae_sedis){
    DF$POWO_complient <-  1
    
    paftol_incertae_sedis_index <- grep("Incertae Sedis", DF$Order)
    taxonomic_fields <- c("Order", "Family", "Genus", "Species")
    
    for (i in paftol_incertae_sedis_index){
      DF$POWO_complient[i] <- 0
       notes <- DF$Taxonomical_Notes[i]
       incertae_sedis_loc <- gregexpr("\\(o\\:.*\\)", notes)
       incertae_sedis_str <- substring(notes,
       									incertae_sedis_loc[[1]] + 1,
       									attr(incertae_sedis_loc[[1]],
       									"match.length") - 1)
       incertae_sedis_cats <- str_split_fixed(unlist(str_split(incertae_sedis_str, ",")),":",2)
       
       DF[i, taxonomic_fields] <- ""
       for (j in 1:dim(incertae_sedis_cats)[1]){
         rank_j <- taxonomic_fields[j]
         name_j <- incertae_sedis_cats[j, 2]
         levels_j <- levels(DF[,rank_j])
         if(!name_j %in% levels_j){
           levels(DF[,taxonomic_fields[j]]) <- c(levels_j, name_j)
         }
         DF[i, rank_j] <- name_j
       }
       
    }
  }
  
  DF$Enriched <- FALSE
  DF$Enriched[DF$SequencingStrategy=="Enriched"] <- TRUE
  DF$Enriched[is.na(DF$Enriched)] <- FALSE
  
  #Create a field for scientific name
  DF$label <- paste0("s", sprintf("%05d",DF$idSequencing))
  DF$Sci_name <- paste0(DF$Genus, "_", DF$Species)
  DF_dims <- dim(DF)
  message(paste0("The initial dataset contains ",
  				DF_dims[1],
  				" lines and ",
  				DF_dims[2],
  				" columns.\n"))
  # Subset data based on a list of ids
  if(!is.null(filter_ids)){
    if(!is.vector(filter_ids)){
      warning("The filter ids must be a vector")
      stop()
    }
    filtering_field_index <- match(filtering_field, colnames(DF))
    if(is.na(filtering_field_index)){
      warning("The filtering field is not present in the table")
      stop()
    }
    
    DF <- DF[DF[,filtering_field_index] %in% filter_ids,]
  }
  
  # Subset data based on a list of columns to output
  if(!is.null(fields)){
    if(!is.vector(fields)){
      warning("The fields element must be a vector")
      stop()
    }
    selected_fields_index <- na.omit(match(fields, colnames(DF)))
    if(length(selected_fields_index) == 0){
      warning("Selected fields don't match column names in the table")
      stop()
    }
    DF <- DF[, selected_fields_index]
  }
  
  
  if(add_clades){
    clades_names <- get_clades_names()
    DF <- cbind(clades_names[match(DF$Family, clades_names$Family), -c(4, 5)], DF)
  }
  
  if(!is.null(recovery_file)){
    recovery <- read.csv(recovery_file, header=T)
    recovery <- recovery[order(-recovery[, 3], -recovery[,2]), ]
    DF$hyb_genes <- recovery[match(DF$label, recovery$sample), 2]
    DF$hyb_genes[is.na(DF$hyb_genes)] <- -1
    DF$hyb_length <- recovery[match(DF$label, recovery$sample), 3]
    DF$hyb_length[is.na(DF$hyb_length)] <- -1
  }

  # Change metadata
  if(!is.null(correct_metadata)){
    new_metadata<- read.csv(correct_metadata, header=T, stringsAsFactors = FALSE)
    for (i in 1:dim(new_metadata)[1]){
      i_line <- new_metadata[i, ]
      if (!(i_line$Target_field %in% colnames(DF) && i_line$Filtering_field %in% colnames(DF) )){
        warning(paste0("Unidentified field in item ", i, ": ", paste(i_line, collapse = "|")))
        next()
      }
      entries_to_change <- DF[, i_line$Filtering_field] == i_line$Filtering_value
      
      DF[, i_line$Target_field] <- add_to_factor(DF[, i_line$Target_field],
      											i_line$Target_value,
      											entries_to_change)
      if(i_line$Target_field == "Sci_name"){
        gen_sp <- unlist(strsplit(i_line$Target_value, " "))
        DF$Genus <- add_to_factor(DF$Genus, gen_sp[1], entries_to_change)
        DF$Species <- add_to_factor(DF$Species, gen_sp[2], entries_to_change)
      }
    }
  }
  
  droplevels(DF)
  DF_dims <- dim(DF)
  message(paste0("The exported dataset contains ",
  				DF_dims[1],
  				" lines and ",
  				DF_dims[2],
  				" columns.\n"))
  save_object(DF, out_file)
}

#### Select best samples (highest recovery) per group

get_best_samples <- function(DF,
                             n = 5,
                             grouping_field = "Family",
                             criterion_field = "hyb_length",
                             out_file = NULL){
    
  DF <- DF[!is.na(DF[,criterion_field]), ]
  
  DF_ordered <- DF[order(DF[ ,criterion_field], decreasing = TRUE), ]
  DF_reduced <- Reduce(rbind, by(DF_ordered, DF_ordered[grouping_field], head, n = n))
  droplevels(DF_reduced)
  save_object(DF_reduced, out_file)
}

#### Produce taxon list ####

taxon_to_tree <- function(DF,
                          drop_blacklisted = TRUE,
                          out_file = NULL){
  
  DF <- DF[!is.na(DF$hyb_genes), ]
  
  if(drop_blacklisted){
    DF <- DF[DF$Blacklisted == 0, ]
  }
  DF <- DF[order(DF$hyb_length, decreasing = TRUE),]
  DF_best_match <- DF[match(unique(DF$idPaftol), DF$idPaftol),]
  
  save_object(DF, out_file)
  
}

#### Produce WCVP Summary ####

# Use a table mapping major clades and orders to families

WCVP_path = "Supporting_files/WCVP_DATA/"
WCVP_sftp = "http://sftp.kew.org/pub/data-repositories/WCVP/"

get_WCVP_summary <- function(WCVP_path = "Supporting_files/WCVP_DATA/",
                             add_clades = FALSE,
                             by_genus = FALSE,
                             update = FALSE,
                             out_file = NULL){
  library(stringr)
  
  # Check if any version is available
  digest_file_suffix <- "_genera.txt"
  digest_files <- dir(path = WCVP_path, pattern = digest_file_suffix)
  latest_digest <- rev(digest_files)[1]

  download <- FALSE
  
  # Flag download if none or not the latest version of WCVP is available
  if(is.na(latest_digest) || update){
    library(RCurl)
    WCVP_releases <- getURL(WCVP_sftp, verbose=FALSE, ftp.use.epsv=TRUE, dirlistonly = TRUE)
    latest_release_file <- rev(getHTMLLinks(WCVP_releases))[1]
    latest_release_version <- sub(".zip", "", latest_release_file)
    latest_release <- paste0(latest_release_version, digest_file_suffix)
    
    if(is.na(latest_digest) || latest_digest != latest_release){
      download <- TRUE
      latest_digest <- latest_release
    } else {
      message("The latest version of WCVP is already in use")
    }
  }
  
  # Download and process latest WCVP version from Kew SFTP
  if (download){
    warning(paste0("Downloading ", latest_release_version))
    latest_release_txt <- paste0(latest_release_version, ".txt")
    orig_wd <- getwd()
    
    setwd(WCVP_path)
    system(paste0("curl ", WCVP_sftp, latest_release_file, " --output ", latest_release_file))
    system(paste0("unzip ", latest_release_file))
    system(paste0("grep '|SPECIES|Accepted|' ",
    				latest_release_txt,
    				" | cut -d '|' -f 2,3 | uniq -c > ",
    				latest_digest))
    system(paste("rm", latest_release_txt, latest_release_file, sep = " "))
    setwd(orig_wd)
  }
  
  wcvp_raw <- read.table(paste0(WCVP_path, "/", latest_digest))
  wcvp_data_species <- data.frame(matrix(str_split_fixed(wcvp_raw[,2], '\\|', 2), ncol = 2)
  								, wcvp_raw[,1])
  colnames(wcvp_data_species) <- c("Family", "Genus", "species_count")
  
  # Summarise by number of genera
  if(by_genus){
    wcvp_data_genus <- data.frame(table(wcvp_data_species$Family))
    colnames(wcvp_data_genus) <- c("Family", "genus_count")
    sp_count <- aggregate(wcvp_data_species$species_count~wcvp_data_species$Family, FUN = sum)
    wcvp_data_genus$species_count <- sp_count[match(wcvp_data_genus$Family, sp_count[,1]), 2]
    wcvp_data <- wcvp_data_genus
  } else {
    wcvp_data <- wcvp_data_species
  }
  
  # Append major clades and orders to the output table
  if(add_clades){
    clades_names <- get_clades_names()
    wcvp_data <- cbind(clades_names[match(wcvp_data$Family, clades_names$Family),-5],
    					 wcvp_data)
  }
  
  # Write output to file or R object
  
  save_object(wcvp_data, out_file)
}

create_seq_list <- function(DF,
							fieltering_field,
							filtering_values,
							full_path = NULL,
							out_file = NULL){
							
  DF_subset <- DF[DF[,fieltering_field] %in% filtering_values,]
  DF_reduced <- get_best_samples(DF_subset, n = 1, grouping_field = "idPaftol")
  labels <- DF_reduced[, "label"]
  if(!is.null(full_path)){
    files_path <- paste0(full_path, "/", substr(labels, 1, 3), "/", labels, ".fasta")
  } else {
    files_path <- paste0(substr(labels, 1, 3), "/", labels, ".fasta")
  }
  save_object(files_path, out_file)
}


####  Get clades names  ####

# Use a csv reference table to map major clades and orders to families
get_clades_names <- function(path = "Supporting_files/clades_orders_families.csv"){
  clades <- read.csv(path)
  clades$Angiosperms <- as.numeric(!clades$Clade_1_rank %in% c("Ferns", "Gymnosperms"))
  clades
}

clades_names <- function(df, path = "Supporting_files/clades_orders_families.csv"){
  clades <- read.csv(path)
  clades$Angiosperms <- as.numeric(!clades$Clade_1_rank %in% c("Ferns", "Gymnosperms"))
  df_appendend <- cbind(clades[match(df$Family, clades$Family),], df)
  df_appendend
}

concat_clades_names <- function(DF){
  
  apply(DF[,c("Clade_1_rank", "Clade_2_rank", "Clade_3_rank")], 1, 
  		function(x) paste(x[x != ""], collapse = "_"))
}

#### Parse node support ####

parse_support <- function(tree, scf = FALSE){
  
  if(scf == FALSE){
    parsed_labels <- lapply(tree$node.label, parse_astral_node_label)
  } else {
    parsed_labels <- lapply(tree$node.label, parse_astral_node_label_scf)
  }
  
  suppressWarnings(parsed_labels_df <- do.call(rbind.data.frame, parsed_labels))
  for (item in 1:dim(parsed_labels_df)[2]){
    tree[[colnames(parsed_labels_df)[item]]] <- as.numeric(parsed_labels_df[,item])
  }
  tree$node <- tree$edge[1,1]:(tree$edge[1,1]+tree$Nnode-1)
  
  return(tree)
}

parse_astral_node_label <- function(node_label, sep = ";"){
  # Consider changing to strsplit_fixed
  node_label <- gsub("'\\[|\\]'", "", node_label)
  items <- unlist(strsplit(unlist(strsplit(node_label, sep)), "="))
  n <- length(items)
  if(n > 1){
    if(n > 4 && items[n-3] == items[n-1]){
      items[n-3] <- "gcf"
    }
    new_labels <- matrix(as.numeric(items[seq(2,n,2)]), 1, n/2, 
    					dimnames = list("", items[seq(1, n-1, 2)]))
    #colnames(new_labels) <- items[seq(1,n-1,2)]
  } else {
    new_labels = ""
  }
  return(new_labels)
}

#### Descendants as names ####

Descendants_as_names <- function(tree, node){
  return(tree$tip.label[unlist(Descendants(tree, node, "tips"))])
}


#### Read Multiple tree into List ####

pyhlos_to_multiphylo <- function(trees_path,
                                 pattern = "",
                                 drop_tips = NULL,
                                 add_names = TRUE,
                                 out_file = NULL){
  library(ape)
  trees_fns <- dir(trees_path, pattern = pattern)
  trees_list <- lapply(paste0(trees_path, "/", trees_fns), FUN = read.tree)
  if(!is.null(drop_tips))
    trees_list <- lapply(trees_list, drop.tip, tip = drop_tips)
  if(add_names){
    names(trees_list) <- trees_fns
  }
  class(trees_list) <- "multiPhylo"
  save_object(trees_list, out_file)
}


#### Sampling on tree ####

get_tree_sampling <- function(tree_fn,
                              by_genus = TRUE,
                              fix_max_ratio = TRUE,
                              out_file = NULL){
  
  library(ape)
  tree <- read.tree(tree_fn)
  
  wcvp_data <- get_WCVP_summary(add_clades = TRUE,
                                   update = TRUE,
                                   by_genus = TRUE)
  
  if(by_genus){
    paftol_data <- unique(export_metadata(filter_ids = as.numeric(tree$tip.label),
                                          fields = c("Order", "Family", "Genus")))
    tree_genus_count <- data.frame(table(paftol_data$Family))
    wcvp_data$tree_genus_count <- tree_genus_count[match(wcvp_data$Family, 
    												tree_genus_count[,1]), 2]
    wcvp_data$tree_genus_count[is.na(wcvp_data$tree_genus_count)] <- 0
    wcvp_data$genus_sampling_ratio <- wcvp_data$tree_genus_count / wcvp_data$genus_count
    if(fix_max_ratio){
      wcvp_data$genus_sampling_ratio[wcvp_data$genus_sampling_ratio > 1] <- 1
    }
    
  } else {
    paftol_data <- unique(export_metadata(filter_ids = as.numeric(tree$tip.label),
                                          fields = c("Order", "Family", "Genus", "Species")))
    wcvp_data$tree_species_count <- tree_genus_count[match(wcvp_data$Family,
    												tree_genus_count[,1]), 2]
    wcvp_data$tree_species_count[is.na(wcvp_data$tree_species_count)] <- 0
    wcvp_data$species_sampling_ratio <- wcvp_data$tree_species_count / wcvp_data$genus_count
    if(fix_max_ratio){
      wcvp_data$species_sampling_ratio[wcvp_data$species_sampling_ratio > 1] <- 1
    }
  }
  
  save_object(wcvp_data, out_file)
}

#### Node age ####

node.age <- function(tree,
                     node){
  node.depths <- node.depth.edgelength(tree)
  node.depths[1] - node.depths[node]
}



#### Ladderise tree ####

ladderise_tree <- function(tree_fn,
                           right = TRUE,
                          out_file = NULL){
  library(ape)
  tree <- read.tree(tree_fn)
  tree_ladderised <- ladderize(tree, right = right)
  
  save_object(tree_ladderised, out_file)
}

#### Rename tree tips ####

rename_tree_tips <- function(tree, DF, in_field, out_field, out_file = NULL){
  if(is.character(tree) && file.exists(tree)){
    tree_obj <- read.tree(tree)
  } else if (class(tree) != "phylo"){
    stop("Tree must be a phylo class object or the path for a newick tree file")
  } else {
    tree_obj <- tree
  }
  new_tip_label <- DF[match(tree_obj$tip.label, DF[, in_field]), out_field]
  if(any(is.na(new_tip_label))){
    is_na_count <- table(is.na(new_tip_label))["TRUE"]
    missing_tips <- tree_obj$tip.label[is.na(new_tip_label)]
    warning(paste0("The following ", 
    		is_na_count, " tip(s) didn't match any entry in the reference dataframe: ",
    		paste(missing_tips, collapse = ", ")))
  }
  tree_obj$tip.label <- new_tip_label
  save_object(tree_obj, out_file)
}

#### Save tip label ####

save_ordered_tip_label <- function(tree_path,
                             out_file = NULL){
  library(ape)
  tree <- read.tree(tree_path)
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  tips_order <- tree$edge[is_tip, 2]
  ordered_tips <- tree$tip.label[tips_order]
  print(class(ordered_tips))
  save_object(ordered_tips, out_file)
}


#### Append Sorta Date Result file ####

append_sorta_data_results <- function(gg_file,
                                      species_tree_path,
                                      gene_trees_path,
                                      gene_tree_pattern = "",
                                      out_file = NULL,
                                      out_sep = ","){
  library(ape)
  
  gg <- read.delim(gg_file, sep=" ", header = TRUE)
  species_tree <- read.tree(species_tree_path)
  gene_trees <- pyhlos_to_multiphylo(gene_trees_path,
  									pattern = gene_tree_pattern, 
  									add_names = TRUE)
  
  Nnode_index <- match("Nnode", names(gene_trees[[1]]))
  gene_trees_Nnode <- unlist(lapply(gene_trees, "[[", Nnode_index))
  
  gg$Nnode <- gene_trees_Nnode[as.character(gg$name)]
  gg$Nnode_rel <- gg$Nnode/species_tree$Nnode
  
  gg$bipartition_rel <- round(gg$bipartition / gg$Nnode_rel, 4)
  
  save_object(gg, out_file, out_sep)
}

simplify_tree <- function(tree,
                          DF,
                          tip_field = "label",
                          summary_field = "FAMILY_PRINT",
                          outliers = NULL,
                          exclude_outliers = TRUE,
                          outliers_label = "short_label",
                          out_file = NULL) {
  
  if(!summary_field %in% colnames(DF)){
    stop("Summary field absent in data frame")
  }
  
  if(length(outliers) > 0){
    if(exclude_outliers){
      excess_tips <- intersect(tree$tip.label, outliers)
      tree <- drop.tip(tree, excess_tips)
    } else {
      tree$tip.label[match(outliers, tree$tip.label)] <- DF[match(outliers, DF[, tip_field]),
      													outliers_label]
    }
  }
  
  tree_DF <- droplevels(DF[DF[, tip_field] %in% tree$tip.label, ])
  no_outiliers_DF <- droplevels(tree_DF[!tree_DF[, tip_field] %in% outliers, ])
  summary_levels <- unique(no_outiliers_DF[, summary_field])
  
  for (level in summary_levels){
    level_tips <- no_outiliers_DF[no_outiliers_DF[, summary_field] == level, tip_field]
    #print(list(level,level_tips))
    if(length(level_tips) > 1){
      tree <- drop.tip(tree, level_tips[-1])
    }
    tree$tip.label[tree$tip.label == level_tips[1]] <- level
  }
  save_object(tree, out_file)
}

#### Prune tree based on data frame

prune_tree <- function(tree_path,
                       taxon_list,
                       new_branch_lengths = FALSE, 
                       out_file = NULL){
  library(ape)
  tree <- read.tree(tree_path)
  
  if (file.exists(as.character(taxon_list))){
    taxon_set <- as.character(read.table(taxon_list))
  } else {
    taxon_set <- as.character(taxon_list)
  }
  
  taxa_to_drop <- tree$tip.label[!tree$tip.label %in% taxon_set]
  pruned_tree <- drop.tip(tree, taxa_to_drop)
  
  if(new_branch_lengths>0){
    pruned_tree$edge.length <- rep(new_branch_lengths, length(pruned_tree$edge.length))
  }
  
  save_object(pruned_tree, out_file = out_file)
}


#### Prepare branch length estimation ####

prepare_branch_length_estimation <- function(gg_file,
                                             gene_thresholds = c(25, 50, 75, 100, 200, 353),
                                             species_tree_path,
                                             gene_trees_path,
                                             gene_tree_pattern,
                                             gene_aligns_path,
                                             gene_aligns_pattern,
                                             sort_by_relative_bp = FALSE,
                                             constraint_to_smaller = FALSE){
  library(ape)
  
  gg <- read.delim(gg_file, sep=" ", header = TRUE)
  
  if(sort_by_relative_bp){
    gg <- gg[order(gg$bipartition_rel, decreasing = TRUE), ]
  }
  
  species_tree <- read.tree(species_tree_path)
  gene_trees <- pyhlos_to_multiphylo(gene_trees_path, pattern = gene_tree_pattern, add_names = TRUE)
  gene_trees_tips_index <- match("tip.label", names(gene_trees[[1]]))
  max_genes <- dim(gg)[1]
  gene_thresholds[gene_thresholds > max_genes] <- max_genes
  gene_thresholds <- sort(unique(gene_thresholds))
  
  gene_trees_th_tips_list <- list()
  
  for (i in 1:length(gene_thresholds)){
    th <- gene_thresholds[i]
    sub_name <- paste0("top_", th)
    gg_th <- gg[1:th, ]
    gene_trees_th <- gene_trees[match(gg_th$name, names(gene_trees))]
    th_gene_names <- gsub(gene_tree_pattern, "", as.character(names(gene_trees_th)))
    gene_trees_th_tips_list[[i]] <- unique(unlist(lapply(gene_trees_th, "[[", gene_trees_tips_index)))
    
    if(sort_by_relative_bp){
      sub_name <- paste0(sub_name, "_rel")
    }
    
    if(constraint_to_smaller){
      gene_trees_th_tips <- gene_trees_th_tips_list[[1]]
      sub_name <- paste0(sub_name, "_small")
    } else {
      gene_trees_th_tips <- gene_trees_th_tips_list[[i]]
    }
    th_missing_tips <- species_tree$tips[!species_tree$tips %in% gene_trees_th_tips]
    species_tree_th <- drop.tip(species_tree, th_missing_tips)
    species_tree_th$edge.length <- rep(0.1, length(species_tree_th$edge.length))
    
    if(dir.exists(sub_name)){
      unlink(sub_name, TRUE)
    }
    
    dir.create(paste0(sub_name, "/genes"), recursive = TRUE)
    write.tree(species_tree_th, paste0(sub_name,"/", sub_name, "_", "constraint.nwk"))
    
    for (gene in th_gene_names){
      print(paste0(gene_aligns_path, "/", gene, gene_aligns_pattern))
      align_in <- read.FASTA(paste0(gene_aligns_path, "/", gene, gene_aligns_pattern))
      align_out <- align_in[names(align_in) %in% gene_trees_th_tips]
      write.FASTA(align_out, paste0(sub_name, "/genes/", gene, ".fasta"))
    }
  }
}  
  
#### Prepare branch length estimation from alignments ####

prepare_branch_length_estimation_from_align <- function(gg_file,
                                             gene_thresholds = c(25, 50, 75, 100, 200, 353),
                                             species_tree_path,
                                             gene_tree_pattern,
                                             gene_aligns_path,
                                             gene_aligns_pattern,
                                             sort_by_relative_bp = FALSE,
                                             constraint_to_smaller = FALSE,
                                             remove_node_labels = TRUE){
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
    dir.create(paste0(sub_name, "/genes"), recursive = TRUE)
    write.tree(species_tree_th, paste0(sub_name,"/", sub_name, "_", "constraint.nwk"))
    
    if(constraint_to_smaller){
      for (j in 1:th){
        align_in <- aligns_list_th[[j]]
        align_out <- align_in[names(align_in) %in% align_th_tips]
        write.FASTA(align_out, paste0(sub_name, "/genes/", genes[j], ".fasta"))
      }
    } else {
      for (j in 1:th){
        write.FASTA(aligns_list_th[[j]], paste0(sub_name, "/genes/", genes[j], ".fasta"))
      }
    }
  message(paste0("Threshold ", th, " prepared."))
  }
}  



#### Get unique taxa ####

get_unique_taxa <- function(df,
                            values,
                            query_field = "idSequencing",
                            return_field){
  
  field_index <- match(query_field, colnames(df))
  taxa <- df[df[, field_index] %in% values, return_field]
  as.character(unique(taxa))
}

