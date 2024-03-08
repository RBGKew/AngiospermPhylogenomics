# LTT heatmap plots 
#
# LTT heatmap plots is a R script for diggesting LTT plots to plot LLT as rasters.
# 
# This script uses a time calibrated tree and  ape::ltt.plot.coords to capture the LTT 
# data to convert into a heatmap, with the colour pattern highlighting periods of
# higher diversification.
# 
# 
# #### Version
# 
# - version = 1.1
# - date: 2024.1.14
# - Currently implementation uses the Lagged Differences of the logaritmic cummulative 
#   sum of histogram count
#
# #### Author
# 
# - Alexandre R. Zuntini
# 

#### Standards #####

LTT_na_value <- "white"
ordinal_LLT_threshold <- 3
age_interval <- 5

##### Packages #####

library(ape)
library(phangorn)
library(phytools)
library(stringr)
library(graphics)
library(ggplot2)
library(reshape2)
library(ggtree)
library(cowplot)
library(gridExtra)
library(readxl)
library(deeptime)
library(caTools)
library(gginnards)
library(rvcheck)
library(aplot)

##### Functions #####

# Extracts order based extra field ($angio_order) in tree object

extract_order_tree <- function(tree, order){
  order_MRCA <- getMRCA(tree, na.omit(tree$tip.label[tree$angio_order == order]))
  order_tree <- extract.clade(tree, node = order_MRCA)
  order_tree$stem_age <- -node.age(tree, Ancestors(tree, order_MRCA, "parent"))
  order_tree
}

#_DEV_#  Add threshold
# Digests LTT plot and create raster overlayed by LTT plot
# Colours must be native to be plotted. NA or something else to skip that element
LTT_heatmap <- function(tree,
                        breaks = 100,
                        smooth = "none",
                        plot_name = NULL,
                        ltt_max_age = 0,
                        cumsum_log_col = "black",
                        ltt_col = "white",
                        LTT_na_value = "gray80",
                        skip_first = TRUE){
  
  ltt_df <- as.data.frame(ltt.plot.coords(tree))
  ltt_df$N_log <- log(ltt_df$N)
  ltt_hist <- hist(ltt_df$time, breaks = breaks, plot = FALSE)
  first_occurrence <- match(TRUE, ltt_hist$counts > 0) # Captures the first diversification as the first occurrence of a non-zero count in the histogram
  
  hist_df <- data.frame(mids = ltt_hist$mids, counts = ltt_hist$counts)
  hist_df$cumsum <- cumsum(hist_df$counts)
  first_cladogenesis <- match(TRUE, hist_df$cumsum > 1)
  hist_df$cumsum_log <- log(hist_df$cumsum)
  hist_df$cumsum_log[1:(first_cladogenesis-1)] <- 0
  hist_df$cumsum_log_diff <-  c(NA, diff(hist_df$cumsum_log))
  if(skip_first){
    hist_df$cumsum_log_diff[first_cladogenesis] <- hist_df$cumsum_log_diff[first_cladogenesis] - log(2)
  }
  diffs <- diff(log(cumsum(hist_df$counts)))
  diffs[1:(first_occurrence - 1)] <- 0
  # Creates the column "diff_log" depending on the smoothing parameter
  hist_df$cumsum_log_diff_mean <- c(NA, runmean(diffs, 3))
  hist_df$cumsum_log_diff_median <- c(NA, smooth(diffs, "S"))
  # Create min and max values for exponential line
  mids_min_max <- summary(hist_df$mids[first_occurrence:length(hist_df$mids)])[c(1, 6)]
  cumsum_log_min_max <- c(0, max(hist_df$cumsum_log))

  if(!is.null(plot_name)){
    p1 <- ggplot(hist_df, aes(x = mids, fill = cumsum_log_diff))
    p1 <- p1 + geom_tile(aes(y = cumsum_log_min_max[2] / 2), alpha = 0.7, height = max(ltt_df$N_log)*2)
    if (cumsum_log_col %in% colors()) {
       p1 <- p1 + geom_point(aes(x = mids, y = cumsum_log),
       						inherit.aes = FALSE,
       						colour = cumsum_log_col)
    }
    if (ltt_col %in% colors()) {
      p1 <- p1 + geom_line(data = ltt_df, aes(x = time, y = N_log),
      						inherit.aes = FALSE,
      						colour = ltt_col)
    }
    p1 <- p1 + scale_fill_viridis_c(na.value = LTT_na_value)
    p1 <- p1 + theme_minimal()
    ggsave(plot_name, p1)
  } else {
    return(hist_df)
  }
}

heatmap_to_bar <- function(df, reverse = FALSE, na.rm = TRUE, max_value = NULL){
  if(reverse){
    df <- df[nrow(df):1, ]
  }
  variable_limits <- c(min(df$variable), 0)
  if(na.rm){
    df <- df[!is.na(df$value), ]
  }
  p <- ggplot(df, aes(x = node, y = variable, fill = value))+
    geom_col(position = position_dodge()) + 
    coord_flip() + 
    ylim(variable_limits) +
    theme_inset()
  if(!is.null(max_value)){
    p <- p + scale_fill_viridis_c(values = c(0, max_value), na.value = LTT_na_value)
  } else {
    p <- p + scale_fill_viridis_c(na.value = LTT_na_value)
  }
  p
}

inset_mod <- function (tree_view, insets, width, height, hjust = 0, vjust = 0, 
          x = "node", reverse_x = FALSE, reverse_y = FALSE) {
  if (height < 0 || height > 1) 
    stop("height should be in range of (0,1)")
  df <- tree_view$data[as.numeric(names(insets)), ]
  x <- match.arg(x, c("node", "branch", "edge"))
  if (x == "node") {
    xx <- df$x
  }
  else {
    xx <- df$branch
  }
  yy <- df$y
  xx <- xx - hjust
  yy <- yy - vjust
  if (reverse_x) 
    xx <- -xx
  if (reverse_y) 
    yy <- -yy
  width <- width * diff(range(tree_view$data$x, na.rm = TRUE))
  height <- height * diff(range(tree_view$data$y, na.rm = TRUE))
  geom_subview <- get_fun_from_pkg("ggimage", "geom_subview")
  tree_view + geom_subview(subview = insets, width = width, 
                           height = height, x = xx, y = yy)
}


##### Data #####

setwd("~/Desktop/LTT/")

source("AngiospermPhylogenomics/Scripts/functions.R")

metadata <- data.frame(read_excel("Figures_tables_files/Supplementary_tables.xlsx", skip = 2))
angio_tree <- read.tree("Trees/4_calibrated_trees/4_young_tree_smoothing_10_pruned_for_diversification_analyses.tre")
angio_tree <- rename_tree_tips(angio_tree, metadata, "Tip_in_tree", "label")
ordinal_outliers <- metadata$label[grep("O", metadata$OUTLIER)]
angio_tree <- drop.tip(angio_tree, ordinal_outliers)
simplified_angio_tree <- ladderize(simplify_tree(angio_tree, metadata, summary_field = "Order"))

tree_metadata <- droplevels(metadata[metadata$label %in% angio_tree$tip.label, ])
unique_genera <- unique(tree_metadata[, c("Order", "Genus")])
order_sampling_df <- data.frame(table(unique_genera$Order))
colnames(order_sampling_df) <- c("label", "sampling")

genera_per_family <- read.csv("Supporting_files/genus_count_per_family.csv")
genera_per_family$Order <- tree_metadata$Order[match(genera_per_family$Family,
													tree_metadata$Family)]
total_genera <- aggregate(genera_per_family$genus_count,
							by = list(genera_per_family$Order),
							FUN = sum)
order_sampling_df$total <- total_genera[match(order_sampling_df$label, total_genera[,1]), 2]
order_sampling_df$new_label <- paste0(order_sampling_df$label,
										" (",
										order_sampling_df$sampling,
										"/",
										order_sampling_df$total,
										")")

# Adds order to tree object. Tree cannot be changed (e.g., rooting or prunning) after this,
# as this field is not modified properly. If changed, the following command must be re-run.
angio_tree$angio_order <- metadata$Order[match(angio_tree$tip.label, metadata$label)]

# Creates series of breaks to be used across Angiosperms orders
angio_breaks <- seq(-160, 0, age_interval)

##### Tree #####

orders <- unique(angio_tree$angio_order)

orders_df <- data.frame(order = orders)
orders_diff <- data.frame(NULL)
retained_orders <- NULL

for (i in 1:dim(orders_df)[1]){
  order <- orders_df$order[i]
  print(order)
  tips <- tree_metadata$label[tree_metadata$Order == order]
  orders_df$ntaxa[i] <- length(tips)
  orders_df$monophyly[i] <- is.monophyletic(angio_tree, tips)
  
  if(length(tips) > 1){
    order_tree <- extract_order_tree(angio_tree, order)
    orders_df$crown_age[i] <- -max(node.age(order_tree))
  } else {
    orders_df$crown_age[i] <- NA
  }
  if(length(tips) >= ordinal_LLT_threshold){
  retained_orders <- c(retained_orders, order)
  diff_log <- LTT_heatmap(order_tree, 
             breaks = angio_breaks,
             smooth = "mean")$cumsum_log_diff
  orders_diff <- rbind(orders_diff, diff_log)
  }
}

for (i in 1:dim(orders_diff)[1]){
  first_occurence <- match(TRUE, orders_diff[i, ] > 0)
  orders_diff[i, 1:(first_occurence-1)] <- NA

}

orders_diff <- cbind(retained_orders, orders_diff)
colnames(orders_diff) <- c("node", angio_breaks[-1])
orders_diff_melt <- melt(orders_diff)
#orders_diff_melt$variable <- -as.numeric(orders_diff_melt$variable)
max_orders_diff_value <- max(orders_diff_melt$value, na.rm = TRUE)
orders_diff_melt$variable <- as.numeric(as.character(orders_diff_melt$variable))
orders_diff_melt_split <- split(orders_diff_melt, orders_diff_melt$node)

nodebars <- lapply(orders_diff_melt_split, heatmap_to_bar) #, max_value = max(orders_diff_melt$value, na.rm = TRUE)
names(nodebars) <- as.character(match(names(nodebars), simplified_angio_tree$tip.label))

# Overall angiosperms 
angio_diffs <- data.frame(node = "Angiosperms",
						 variable = angio_breaks[-1], 
						 value = LTT_heatmap(angio_tree, 
                         breaks = angio_breaks,
                         smooth = "mean")$cumsum_log_diff)

tree_xlim = 40
root.position = -154
label_rel_size <- tree_xlim/(tree_xlim - root.position)

# Tree with heatmap
p_tree <- ggtree(simplified_angio_tree, root.position = root.position)
p_tree <- p_tree  %<+% order_sampling_df + geom_tiplab(mapping = aes(label = new_label), size = 2) +
  xlim_tree(tree_xlim) +
  theme_tree(plot.margin = margin(0, 0, 0, 0))

orders_df_clean <- na.omit(orders_df)
colnames(orders_df_clean)[1] <- "label"

p_tree_nodebars <- inset_mod(p_tree,
							nodebars, width = 1.115,
							height = 0.022,
							hjust = 78.1,
							vjust = 0.2)
p_tree_nodebars <- p_tree_nodebars %<+% orders_df_clean + 
                    geom_text(aes(x = crown_age), 
                    label = "|",
                    colour = "black",
                    nudge_y = 0,
                    size = 3,
                    fontface = 2)

p_angio <- heatmap_to_bar(angio_diffs) + 
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))

p_blank <- ggplot() + theme_nothing()

p_geo <-  ggplot() +
  scale_x_continuous(limits = c(-160,0)) +
  theme_inset() +
  theme(axis.text.x = element_text(colour = "black", size = 5),
  						axis.ticks.x = element_line(colour = "black", linewidth = 0.5))+
  coord_geo(pos = as.list(rep("bottom", 2)),
            dat = list("periods", "eras"),
            fill = "white",
            height = list(unit(0.5, "lines"),
            unit(0.5, "line")),
            rot = list(0,  0),
            size = list(1.5, 1.5), 
            abbrv = FALSE,
            neg = TRUE,
            center_end_labels = TRUE)

p_angio_line <- align_plots(p_angio, p_blank, axis = "lr")
p_angio_grid <- grid.arrange(p_angio_line[[1]], p_angio_line[[2]],
							ncol = 2, nrow = 1,
							heights = 1, widths = c(1, 0.24))

p_geo_line <- align_plots(p_geo, p_blank, axis = "lr")
p_geo_grid <- grid.arrange(p_geo_line[[1]], p_geo_line[[2]], 
							ncol = 2, nrow = 1,
							heights = 1, widths = c(1, 0.3))

plots <- align_plots(p_tree_nodebars, p_angio_grid, p_geo_grid, axis = "lr")
gg_final <- grid.arrange(plots[[1]], plots[[2]], plots[[3]],
								ncol = 1, nrow = 3,
								heights = c(6, 0.3, 0.2), widths = 8)
                    
filename <- "Figure_2.pdf"

ggsave(filename, gg_final, width = 6, height = 7.5)
system2("open", filename)

