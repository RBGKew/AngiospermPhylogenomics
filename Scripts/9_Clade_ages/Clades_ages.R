#!/usr/bin/env Rscript
#
# Clades_ages.R
#
# Clades_ages.R retrieves the ages of major clades, orders and families and compares with
# Ramírez-Barahona et al. 2020. Nature Ecology and Evolution 4: 1232-1238.
# 
# 
# #### Version
# 
# - version = 1.0
# - date: 2023.5.12
#
#
# #### Author
# 
# - Alexandre R. Zuntini
# - a.zuntini@kew.org
# 

library(ape)
library(phangorn)
library(readxl)
library(ggplot2)
library(cowplot)
library(gridExtra)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("functions.R")

theme_set(theme_get() + theme(text = element_text(family = 'sans')))

node_age <- function(tree,
                     node,
                     node_depths = NULL,
                     include_stem = FALSE){
  if(is.null(node_depths)){
    node_depths <- node.depth.edgelength(tree)
  }
  crown_node_age <- node_depths[1] - node_depths[node]
  if(include_stem){
    stem_node <- tree$edge[tree$edge[,2] == node, 1]
    if(length(stem_node) == 0){
      stem_node_age <- NA
    } else {
      stem_node_age <- node_depths[1] - node_depths[stem_node]
    }
    c(crown_node_age, stem_node_age)
  } else {
    c(crown_node_age, NA)
  }
}


young_tree <- read.tree("Trees/4_young_tree_smoothing_10_pruned_for_diversification_analyses.tre.tre")
old_tree <- read.tree("Trees/4_old_tree_smoothing_10_pruned_for_diversification_analyses.tre")
metadata <- data.frame(read_excel("Figures_tables_files/Supplementary_Tables.xlsx"))


rb_data_CC_complete <- data.frame(read_excel("Ramirez-Barahona_2020_Data4_Ages.xlsx",
                                             sheet = "CC_complete_Ages",
                                             col_types =c("text",
                                             			rep("numeric", 9),
                                             			rep("text", 3))))
colnames(rb_data_CC_complete)[1] <- "level"
rb_data_UC_complete <- data.frame(read_excel("Ramirez-Barahona_2020_Data4_Ages.xlsx",
                                             sheet = "UC_complete_ages",
                                             col_types =c("text",
                                             			rep("numeric", 9),
                                             			rep("text", 3))))

colnames(rb_data_UC_complete)[1] <- "level"

young_tree <- rename_tree_tips(young_tree,
                            metadata,
                            in_field = "Tip_in_tree", 
                            out_field = "label")
old_tree <- rename_tree_tips(old_tree, 
                               metadata,
                               in_field = "Tip_in_tree", 
                               out_field = "label")

young_tree_no_outliers <- drop.tip(young_tree,
							tip = metadata$label[!is_blank(metadata$OUTLIER)])
old_tree_no_outliers <- drop.tip(old_tree,
							tip = metadata$label[!is_blank(metadata$OUTLIER)])
tree_metadata <- metadata[metadata$label %in% young_tree_no_outliers$tip.label, ]

ranks <- c("Clade_1_rank", "Clade_2_rank", "Clade_3_rank", "Order", "Family")

age_df <- data.frame(NULL)
young_node.depths <- node.depth.edgelength(young_tree_no_outliers)
old_node.depths <- node.depth.edgelength(old_tree_no_outliers)

for (rank in ranks){
  print(rank)
  rank_df <- data.frame(rank = rank,
                        level = na.omit(unique(tree_metadata[, rank])))
  levels_tips <- lapply(rank_df$level,
  				FUN = function(x) na.omit(tree_metadata$label[tree_metadata[, rank] == x]))
  rank_df$number_tips <- unlist(lapply(levels_tips, FUN = length))
  levels_MRCA <- lapply(levels_tips, FUN = getMRCA, phy = young_tree_no_outliers)
  for (i in 1:length(levels_MRCA)){
    if (is.null(levels_MRCA[[i]])){
      levels_MRCA[[i]] <- match(levels_tips[[i]], young_tree_no_outliers$tip.label)
    }
  }
  rank_df$MRCA <- unlist(levels_MRCA)
  rank_df$tips_found <- unlist(lapply(levels_MRCA,
  				FUN = function(x) length(Descendants(young_tree_no_outliers,
  													node = x ,
  													type = "tips"))))
  young_tree_ages <- as.data.frame(do.call(rbind,
                            lapply(levels_MRCA,
                            FUN = node_age, 
                            tree = young_tree_no_outliers,
                            node_depths = young_node.depths,
                            include_stem = TRUE)))
  old_tree_ages <- as.data.frame(do.call(rbind,
                            lapply(levels_MRCA,
                            FUN = node_age, 
                            tree = old_tree_no_outliers,
                            node_depths = old_node.depths,
                            include_stem = TRUE)))
  rank_df <- cbind(rank_df, young_tree_ages, old_tree_ages)
  #rank_df[rank_df$number_tips == 1, c(6,8)] <- NA
  age_df <- rbind(age_df, rank_df)
}
colnames(age_df)[6:9] <- c("Young_tree_crown_node",
							"Young_tree_stem_node",
							"Old_tree_crown_node",
							"Old_tree_stem_node")
write.table(age_df, "clade_ages.txt", sep = "\t")

young_age_combined <- merge(age_df, rb_data_CC_complete,  all.x = TRUE)
old_age_combined <- merge(age_df, rb_data_UC_complete,  all.x = TRUE)

young_age_combined$diff_stem_node <- young_age_combined$Young_tree_stem_node - 
										young_age_combined$Stem_treePL
old_age_combined$diff_stem_node <- old_age_combined$Old_tree_stem_node -
										old_age_combined$Stem_treePL

young_age_orders_families <- droplevels(young_age_combined[young_age_combined$rank %in% 
										c("Order", "Family"), ])
old_age_orders_families <- droplevels(old_age_combined[old_age_combined$rank %in% 
										c("Order", "Family"), ])


plot_points <- function(df, x, y, colour, title, xlim = c(40, 160), filename = NULL){
  sub_df <- df[, c(x, y, colour)]
  colnames(sub_df) <- c("x", "y", "colour")
  p <- ggplot(sub_df, aes(x = x, y = y, colour = colour)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    labs(title = title, x = "Our tree", y = "Ramirez-Barahona 2020") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = xlim) +
    scale_colour_manual(name = NULL, values = c("#D55E00", "#56B4E9")) + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(linewidth = 1, colour = "black"),
          text = element_text(size = 16, colour = "black"),
          axis.ticks = element_line(linewidth = 1, colour = "black"))
  if(is.null(filename)){
    p
  } else {
    ggsave(filename, p)
  }
}

plot_diffs <- function(df, x, y, title, filename = NULL){
  sub_df <- df[, c(x, y)]
  colnames(sub_df) <- c("x", "y")
  p <- ggplot(sub_df, aes(x = x, y = y, fill = x)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_violin() +
    geom_jitter() +
    labs(title = title, x = "Taxonomic rank", y = "Age difference\n(Our tree - Ramirez-Barahona)") +
    scale_fill_manual(breaks = c("Order", "Family"), values = c("#D55E00", "#56B4E9")) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(linewidth = 1, colour = "black"),
          text = element_text(size = 16, colour = "black"),
          axis.ticks = element_line(linewidth = 1, colour = "black")) +
    theme(legend.position = "none")
  if(is.null(filename)){
    p
  } else {
    ggsave(filename, p)
  }
}


young_x_y <- plot_points(young_age_orders_families,
						"Young_tree_stem_node",
						"Stem_treePL",
						"rank",
						title = "a)")+ 
  theme(axis.title.x = element_blank()) +
  theme(legend.position = c(0.15, 0.9))
old_x_y<- plot_points(old_age_orders_families,
						"Old_tree_stem_node",
						"Stem_treePL",
						"rank",
						xlim = c(40, 250),
						title = "c)")  +
  theme(legend.position = "none")

young_diffs <- plot_diffs(young_age_orders_families, "rank", "diff_stem_node", "b)") + 
  theme(axis.title.x = element_blank())
old_diffs <- plot_diffs(old_age_orders_families, "rank", "diff_stem_node", "d)")

young_plots <- align_plots(young_x_y, young_diffs, align = 'h', axis = 'tb')
old_plots <- align_plots(old_x_y, old_diffs, align = 'h', axis = 'tb')

gg_final <- grid.arrange(young_plots[[1]], young_plots[[2]],
							old_plots[[1]], old_plots[[2]],
							ncol = 2, nrow = 2,
							widths = c(0.9, 1), heights = c(0.9, 1))
ggsave("Age_comparsion.pdf", gg_final, width = 10, height = 10)

