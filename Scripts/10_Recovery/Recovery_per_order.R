library(ggplot2)
library(ggtree)
library(scales)
library(reshape2)
library(readxl)
library(ape)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../functions.R")

metadata <- data.frame(read_excel("Figures_tables_files/Supplementary_Tables",
									skip = 2,
									col_types = "text"))
astral_tree <- read.tree("Trees/2_global_tree_original_label.tre")
seqs_lengths <- read.csv("no_gap_length.csv", row.names = 1)
align_lengths <- read.csv("total_length.csv", row.names = 1)

samples_per_order <- data.frame(table(metadata$Order))
colnames(samples_per_order) <- c("label", "samples")
samples_per_order$PRINT <- paste0(samples_per_order$label,
								" (",
								samples_per_order$samples, ")",
								paste0(rep(".", 40), collapse = ""))

outliers <- metadata$label[!is.na(metadata$OUTLIER)]
outliers <- c(outliers, metadata$label[metadata$Clade_1_rank == "Gymnosperms"])
astral_tree$edge.length <- rep(1, length(astral_tree$edge.length))
order_tree <- simplify_tree(astral_tree,
                            metadata,
                            tip_field = "label",
                            summary_field = "Order",
                            outliers = outliers,
                            exclude_outliers = TRUE)
order_tree <- ladderize(root(order_tree, "Amborellales", resolve.root = TRUE))
p_tree <- ggtree(order_tree, ladderize = TRUE)
p_tree <- p_tree %<+% samples_per_order + geom_tiplab(aes(label = PRINT), size = 3) +
  xlim_tree(120)

#### Recovery ####
gene_lengths <- apply(align_lengths, 2, max)
seqs_lengths_rel <- sweep(as.matrix(seqs_lengths), 2, gene_lengths, `/`)
seqs_lengths_rel_NA <- ifelse(seqs_lengths_rel == 0, NA, seqs_lengths_rel)
seqs_lengths_rel_NA <- data.frame(seqs_lengths_rel_NA)
seqs_lengths_rel_NA$Order <- metadata$Order[match(rownames(seqs_lengths_rel_NA),
															metadata$label)]
order_recovery <- aggregate(seqs_lengths_rel_NA,
							by = list(label = seqs_lengths_rel_NA$Order),
							FUN = mean,
							na.rm = TRUE)
rownames(order_recovery) <- order_recovery$label
order_recovery <- order_recovery[, -1]

samples_per_order <- data.frame(table(metadata$Order))
colnames(samples_per_order) <- c("label", "samples")
samples_per_order$PRINT <- paste0(samples_per_order$label,
									" (",
									samples_per_order$samples,
									")",
									paste0(rep(".", 40), collapse = ""))

p_tree_recovery <- gheatmap(p_tree, order_recovery, offset = 17, width = 3.5, 
         colnames=FALSE) +
  scale_y_continuous(expand=c(0, 0.3)) +
  scale_fill_viridis_c(option="D", name="Average\nrecovery", na.value = "white") +
  theme(legend.position = c(0.06, 0.87))
ggsave("Recovery.pdf", p_tree_recovery)

#### Occupancy ####
seqs_occupancy <- data.frame(ifelse(as.matrix(seqs_lengths) > 0, 1, 0))
seqs_occupancy$Order <- metadata$Order[match(rownames(seqs_occupancy), metadata$label)]
order_occupancy <- aggregate(seqs_occupancy,
								by = list(label = seqs_occupancy$Order),
								FUN = mean)
rownames(order_occupancy) <- order_occupancy$label
order_occupancy <- order_occupancy[, -1]

p_tree_occupancy <- gheatmap(p_tree, order_occupancy, offset = 17, width = 3.5, 
                            colnames=FALSE) +
  scale_y_continuous(expand=c(0, 0.3)) +
  scale_fill_viridis_c(option="D", name="Average\noccupancy", na.value = "white") +
  theme(legend.position = c(0.06, 0.87))
ggsave("Occupancy.pdf", p_tree_occupancy)

