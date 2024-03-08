library(ggplot2)
library(coda)
library(phytools)
library(devtools)
library(phylobase)
library(phangorn)
library(data.table)
library(TreeSim)

##################
###OVERALL_INFO###
##################

###overall tree

tree <- read.tree("young_tree_smoothing_10_pruned_for_diversification_analyses.tre") # "young_tree_smoothing_10_pruned_for_diversification_analyses.tre"

###diversification parameters

missing_tree_net_div <- log(360000/2.0) / max(node.depth.edgelength(tree)) # just leave it
missing_tree_ext <- missing_tree_net_div*99     # for 95 is 19, for 99 is 99, for 97.5 is 39
missing_tree_spec <- missing_tree_net_div*100    # for 95 is 20, for 99 is 100, for 97.5 is 40 the difference always has to be one diversification rate unit

################################################
###GET_KEY_INFORMATION_ON_GENERA_AND_FAMILIES###
################################################

###tree tips just as genera or families

genera_only_tips <- vector(mode="numeric", length=0)
families_only_tips <- vector(mode="numeric", length=0)
for (i in 1:length(tree$tip.label)){
genera_only_tips <- append(genera_only_tips, unlist(strsplit(tree$tip.label[[i]], split="_"))[[3]])
families_only_tips <- append(families_only_tips, unlist(strsplit(tree$tip.label[[i]], split="_"))[[2]])
}

###get unique genera and families

genera <- unique(genera_only_tips)
families <- unique(families_only_tips)

#########################################
###GET_INFO_FOR_GENERA_IN_OVERALL_TREE###
#########################################

###get all tip labels for each genus

genera_tips <- vector("list", length(genera))
for (i in 1:length(tree$tip.label)){
genera_tips[[which(genera==unlist(strsplit(tree$tip.label[[i]], split="_"))[[3]])]] <- append(genera_tips[[which(genera==unlist(strsplit(tree$tip.label[[i]], split="_"))[[3]])]], tree$tip.label[[i]])
}

###only one tip per genus - so this equals genera tips

genera_mrca_clades <- vector("list", length(genera))
for (i in 1:length(genera)){
if (length(genera_tips[[i]]) > 1){
genera_mrca_clades[[i]] <- extract.clade(tree, findMRCA(tree, genera_tips[[i]], "node"))$tip.label
} else {
genera_mrca_clades[[i]] <- genera_tips[[i]]
}
}

###again ver simple as only one tip per genus

final_tip_labels_genera <- vector("list", length(genera_tips))
final_node_genera <- vector("list", length(genera_tips))
for (i in 1:length(genera_tips)){
if (length(genera_tips[[i]]) > 0){
if ((unlist(strsplit(genera_tips[[i]], split="_"))[[3]] == "Incertae") == FALSE){ 
if (length(genera_tips[[i]]) > 1){
final_tip_labels_genera[[i]] <- extract.clade(tree, findMRCA(tree, genera_tips[[i]], "node"))$tip.label
final_node_genera[[i]] <- findMRCA(tree, final_tip_labels_genera[[i]], "node")
} else {
final_tip_labels_genera[[i]] <- tree$tip.label[[which(tree$tip.label == genera_tips[[i]])]]
final_node_genera[[i]] <- which(tree$tip.label == final_tip_labels_genera[[i]][[1]])
}
}
}
}

analysis_genera_name <- genera

###UNLIKE WITH FAMILIES BELOW THERE IS NO EXTRA CLEANING PHASE HERE AS ONLY ONE PER GENUS###

###########################################
###GET_INFO_FOR_FAMILIES_IN_OVERALL_TREE###
###########################################

###get all tip labels for each family

families_tips <- vector("list", length(families))
for (i in 1:length(tree$tip.label)){
families_tips[[which(families==unlist(strsplit(tree$tip.label[[i]], split="_"))[[2]])]] <- append(families_tips[[which(families==unlist(strsplit(tree$tip.label[[i]], split="_"))[[2]])]], tree$tip.label[[i]])
}

###get mrca tips of initial family tips

families_mrca_clades <- vector("list", length(families))
for (i in 1:length(families)){
if (length(families_tips[[i]]) > 1){
families_mrca_clades[[i]] <- extract.clade(tree, findMRCA(tree, families_tips[[i]], "node"))$tip.label
} else {
families_mrca_clades[[i]] <- families_tips[[i]]
}
}

###mark families mrca clades that are too big

removal <- vector(mode="numeric", length=0)
for (i in 1:length(families_tips)){
if (((length(families_mrca_clades[[i]]) - length(families_tips[[i]]))/length(families_tips[[i]])) > 0.25){
removal <- append(removal, i)
}
}

###try_to_re_add_missing_by_getting_biggest block_belonging_to_a_taxon

for (j in removal){
blocks <- which(families_only_tips == families[[j]])
block_end_positions <- vector(mode="numeric", length=0)
block_sizes <- vector(mode="numeric", length=0)
for (i in 1:(length(blocks)-1)){
if (blocks[[i]] < blocks[[i+1]]-3){ # new term - allows the insertion of one random tip - allows two inserted tips
block_end_positions <- append(block_end_positions, i)
}
}
block_end_positions <- append(block_end_positions, length(blocks))
block_end_positions <- append(block_end_positions, 0, after=0)
block_sizes <- 0
for (i in 2:length(block_end_positions)){
block_sizes <- append(block_sizes, block_end_positions[[i]] - block_end_positions[[i-1]])
}
if (length(which(block_sizes == max(block_sizes))) == 1){
if (max(block_sizes)>1){
if (max(block_sizes) >= (2/3)*sum(block_sizes)){
tip_label_vector_to_add <- blocks[seq(block_end_positions[[which(block_sizes==max(block_sizes))-1]]+1,block_end_positions[[which(block_sizes==max(block_sizes))]], 1)]  
families_tips[[j]] <- tree$tip.label[tip_label_vector_to_add]
}
}
} else {
families_tips[[j]] <- vector(mode="numeric", length=0) # if unable to redefine taxon set, length set to 0 and removed below
}
}

###GET_FINAL_NODES_AND_TIPS

final_tip_labels_families <- vector("list", length(families_tips))
final_node_families <- vector("list", length(families_tips))
for (i in 1:length(families_tips)){
if (length(families_tips[[i]]) > 0){
if ((unlist(strsplit(families_tips[[i]], split="_"))[[3]] == "Incertae") == FALSE){ 
if (length(families_tips[[i]]) > 1){
final_tip_labels_families[[i]] <- extract.clade(tree, findMRCA(tree, families_tips[[i]], "node"))$tip.label
final_node_families[[i]] <- findMRCA(tree, final_tip_labels_families[[i]], "node")
} else {
final_tip_labels_families[[i]] <- tree$tip.label[[which(tree$tip.label == families_tips[[i]])]]
final_node_families[[i]] <- which(tree$tip.label == final_tip_labels_families[[i]][[1]])
}
}
}
}

###mark 0 length for removal

removal <- vector(mode="numeric", length=0)
for (i in 1:length(final_tip_labels_families)){
if (length(final_tip_labels_families[[i]]) == 0){
removal <- append(removal, i)
} 
}

if (length(removal) > 0){
final_tip_labels_families <- final_tip_labels_families[-removal]
final_node_families <- final_node_families[-removal]
analysis_families_name <- families[-removal]
}

###remove things belonging to another family - i.e. if node defines a grade, or such like

for (a in 1:length(final_tip_labels_families)){
final_tip_labels_families_only_families <- unlist(lapply(lapply(lapply(final_tip_labels_families[[a]], strsplit, "_"), unlist), "[[", 2))
tester <- which(final_tip_labels_families_only_families != analysis_families_name[[a]])
if (length(tester) > 0){
final_tip_labels_families[[a]] <- final_tip_labels_families[[a]][-tester]
print(final_tip_labels_families[[a]][-tester])
}
}

##############################
###DO_THE_SIMULATION##########
##############################

###get other family variables

stem_nodes_families <- vector(mode="numeric", length=0)
stem_times_families <- vector(mode="numeric", length=0)
for (i in 1:length(final_tip_labels_families)){
if (length(final_tip_labels_families[[i]]) > 1){
stem_nodes_families <- append(stem_nodes_families, tree[[1]][,1][[which(tree[[1]][,2] == final_node_families[[i]])]])
stem_times_families <- append(stem_times_families, max(node.depth.edgelength(extract.clade(tree, stem_nodes_families[[length(stem_nodes_families)]]))))
}
if (length(final_tip_labels_families[[i]]) == 1){
stem_nodes_families <- append(stem_nodes_families, tree[[1]][,1][[which(tree[[1]][,2] == final_node_families[[i]])]])
stem_times_families <- append(stem_times_families, max(node.depth.edgelength(extract.clade(tree, stem_nodes_families[[length(stem_nodes_families)]]))))
}
}

###get other genus variables

stem_nodes_genera <- vector(mode="numeric", length=0)
stem_times_genera <- vector(mode="numeric", length=0)
genera_corresponding_families_names <- vector(mode="numeric", length=0)
genera_corresponding_families_stem_time <- vector(mode="numeric", length=0)
for (i in 1:length(final_tip_labels_genera)){
addition <- vector(mode="numeric", length=0)
if (length(final_tip_labels_genera[[i]]) > 1){
if (length(which(analysis_families_name == unlist(strsplit(final_tip_labels_genera[[i]][[1]], split="_"))[[2]])) > 0){ # can create an issue where non mono fams excluded
stem_nodes_genera <- append(stem_nodes_genera, tree[[1]][,1][[which(tree[[1]][,2] == final_node_genera[[i]])]])
stem_times_genera <- append(stem_times_genera, max(node.depth.edgelength(extract.clade(tree, stem_nodes_genera[[length(stem_nodes_genera)]]))))
genera_corresponding_families_names <- append(genera_corresponding_families_names, unlist(strsplit(final_tip_labels_genera[[i]][[1]], split="_"))[[2]])
genera_corresponding_families_stem_time <- append(genera_corresponding_families_stem_time, stem_times_families[[which(analysis_families_name == unlist(strsplit(final_tip_labels_genera[[i]][[1]], split="_"))[[2]])[[1]]]]) # fine to condition on tip here bc one tip in each genus
addition <- append(addition, 1)
} 
}
if (length(final_tip_labels_genera[[i]]) == 1){
if (length(which(analysis_families_name == unlist(strsplit(final_tip_labels_genera[[i]][[1]], split="_"))[[2]])) > 0){ # can create an issue where non mono fams excluded
stem_nodes_genera <- append(stem_nodes_genera, tree[[1]][,1][[which(tree[[1]][,2] == final_node_genera[[i]])]])
stem_times_genera <- append(stem_times_genera, max(node.depth.edgelength(extract.clade(tree, stem_nodes_genera[[length(stem_nodes_genera)]]))))
genera_corresponding_families_names <- append(genera_corresponding_families_names, unlist(strsplit(final_tip_labels_genera[[i]][[1]], split="_"))[[2]])
genera_corresponding_families_stem_time <- append(genera_corresponding_families_stem_time, stem_times_families[[which(analysis_families_name == unlist(strsplit(final_tip_labels_genera[[i]][[1]], split="_"))[[2]])[[1]]]])
addition <- append(addition, 1)
}
}
if (length(addition) == 0){
removal <- append(removal, i)
}
}

final_tip_labels_genera <- final_tip_labels_genera[-removal]
final_node_genera <- final_node_genera[-removal]
analysis_genera_name <- analysis_genera_name[-removal]

###read_in_wcvp_data###
table <- fread("wcvp_v8_mar_2022.txt")
table <- table[which(table[,9][[1]] == "Accepted"),]

###starter variables

intervals <- seq(0, 160, 1)
interval_missing_diversity_rep_storage <- vector("list", length(intervals)) 

###begin

for (z in 1:10){ # number of times to run simulation
interval_missing_diversity_storage <- vector("list", length(intervals))

###DO_ADDITION_OF_MISSING_SPECIES_IN_GENUS###
#############################################

genera_done <- vector(mode="numeric", length=0)
for (a in 1:length(final_node_families)){
number_of_species_done <- vector(mode="numeric", length=0)
if (final_node_families[[a]] > length(tree$tip.label)){ ###if the family is not a tip
families_genera <- vector(mode="numeric", length=0)
for (b in 1:length(final_tip_labels_families[[a]])){ ###number_of_genera_in_families_in_tree
families_genera <- append(families_genera, unlist(strsplit(final_tip_labels_families[[a]][[b]], split="_"))[[3]])
}
families_genera <- unique(families_genera) # get the genera in each family again
removal <- vector(mode="numeric", length=0) # need to remove genera that have inadvertantly been done within another family
if (length(genera_done) > 0){
for (i in 1:length(families_genera)){
if (families_genera[[i]] %in% genera_done){
removal <- append(removal, i)
}
}
}
if (length(removal) > 0){
families_genera <- families_genera[-removal]
}
genera_done <- append(genera_done, families_genera)
if (length(families_genera) > 0){ # if there are still ones left to do
for (b in 1:length(families_genera)){ # get_diversity_of_each_genus_in_family
genus_subset_table <- table[which(table[,3]==families_genera[[b]]),]
genus_diversity <- nrow(genus_subset_table)
message(paste(families_genera[[b]], ", which is in the tree has ", genus_diversity, " species in WCVP", sep="")) 
if (length(which(analysis_genera_name == families_genera[[b]])) == 1){ ###check on whether genus is in the family
if (final_node_genera[[which(analysis_genera_name == families_genera[[b]])]] > length(tree$tip.label)){ ###get_genus_diversity_in_the_tree if not a tip
genus_diversity_in_tree <- length(extract.clade(tree, final_node_genera[[which(analysis_genera_name == families_genera[[b]])]])$tip.label)
} else {
genus_diversity_in_tree <- 1
}
if ((genus_diversity - genus_diversity_in_tree) > 0){ ###simulate_tree_if_conditions_are_met
genus_stem_time <- stem_times_genera[[which(analysis_genera_name == families_genera[[b]])]]
genus_tree <- sim.bd.taxa.age(genus_diversity, 1, missing_tree_spec, missing_tree_ext, frac = 1, genus_stem_time, mrca = FALSE)[[1]]
number_of_species_done <- append(number_of_species_done, length(genus_tree$tip.label))
#message(paste("initial simulated genus tree diversity = ", length(genus_tree$tip.label), sep="")) 
initial_age_storage <- max(node.depth.edgelength(extract.clade(genus_tree, findMRCA(genus_tree, genus_tree$tip.label, "node"))))
genus_tree <- drop.tip(genus_tree, sample(seq(1, length(genus_tree$tip.label), 1), genus_diversity_in_tree)) # drop sufficient tips from simulated tree 
#message(paste("final missing simulated genus tree diversity = ", length(genus_tree$tip.label), sep=""))
if (length(genus_tree$edge.length) == 1){
genus_tree$edge.length[[1]] <- initial_age_storage
}
genus_tree_table <- data.frame(genus_tree[[1]], genus_tree$edge.length, seq(1, nrow(genus_tree[[1]]), 1))
genus_tree_table <- genus_tree_table[order(genus_tree_table[,2]),]
end_times <- c(rep(0, length(genus_tree$tip.label)), branching.times(genus_tree)[-1])
start_times <- end_times + genus_tree_table[,3]
for (c in 1:length(interval_missing_diversity_storage)){
counter <- length(which(end_times <= intervals[[c]] & start_times > intervals[[c]]))
if (counter > 0){
interval_missing_diversity_storage[[c]] <- append(interval_missing_diversity_storage[[c]], counter)
}  
}
for (c in 1:length(interval_missing_diversity_storage)){
if (length(interval_missing_diversity_storage[[c]]) > 0){
interval_missing_diversity_storage[[c]] <- sum(interval_missing_diversity_storage[[c]])
}
}
}
}
}
}
} else {
families_genera <- vector(mode="numeric", length=0)
families_genera <- unlist(strsplit(tree$tip.label[[final_node_families[[a]]]], split="_"))[[3]] ###get_the_genus_of_the_single_tip_in_the_tree
removal <- vector(mode="numeric", length=0) # filter
if (length(genera_done) > 0){ 
for (i in 1:length(families_genera)){
if (families_genera[[i]] %in% genera_done){
removal <- append(removal, i)
}
}
}
if (length(removal) > 0){
families_genera <- families_genera[-removal]
}
genera_done <- append(genera_done, families_genera)
if (length(families_genera) > 0){
for (b in 1:length(families_genera)){ # get_diversity_of_each_genus_in_family
genus_subset_table <- table[which(table[,3]==families_genera[[b]]),]
genus_diversity <- nrow(genus_subset_table)
message(paste(families_genera[[b]], ", which is in the tree has ", genus_diversity, " species in WCVP", sep="")) 
if (length(which(analysis_genera_name == families_genera[[b]])) == 1){ ###check on whether genus is in the family
genus_diversity_in_tree <- 1
if ((genus_diversity - genus_diversity_in_tree) > 0){ ###simulate_tree_if_conditions_are_met
genus_stem_time <- stem_times_genera[[which(analysis_genera_name == families_genera[[b]])]]
genus_tree <- sim.bd.taxa.age(genus_diversity+1, 1, missing_tree_spec, missing_tree_ext, frac = 1, genus_stem_time, mrca = FALSE)[[1]]
number_of_species_done <- append(number_of_species_done, length(genus_tree$tip.label))
#message(paste("initial simulated genus tree diversity = ", length(genus_tree$tip.label), sep="")) 
initial_age_storage <- max(node.depth.edgelength(extract.clade(genus_tree, findMRCA(genus_tree, genus_tree$tip.label, "node"))))
genus_tree <- drop.tip(genus_tree, sample(seq(1, length(genus_tree$tip.label), 1), genus_diversity_in_tree)) #clarify_single_branches - which branching time are you left with 
#message(paste("final missing simulated genus tree diversity = ", length(genus_tree$tip.label), sep=""))
if (length(genus_tree$edge.length) == 1){
genus_tree$edge.length[[1]] <- initial_age_storage
}
genus_tree_table <- data.frame(genus_tree[[1]], genus_tree$edge.length, seq(1, nrow(genus_tree[[1]]), 1))
genus_tree_table <- genus_tree_table[order(genus_tree_table[,2]),]
end_times <- c(rep(0, length(genus_tree$tip.label)), branching.times(genus_tree)[-1])
start_times <- end_times + genus_tree_table[,3]
for (c in 1:length(interval_missing_diversity_storage)){
counter <- length(which(end_times <= intervals[[c]] & start_times > intervals[[c]]))
if (counter > 0){
interval_missing_diversity_storage[[c]] <- append(interval_missing_diversity_storage[[c]], counter)
}  
}
for (c in 1:length(interval_missing_diversity_storage)){
if (length(interval_missing_diversity_storage[[c]]) > 0){
interval_missing_diversity_storage[[c]] <- sum(interval_missing_diversity_storage[[c]])
}
}
}
}
}
}
}
###now do the unsampled genera
if (length(final_node_families[[a]]) > 0){ # it should be by definition of the filtering previously
if (final_node_families[[a]] > length(tree$tip.label)){ ###if it refers to a node rather than a tip
missing_genera <- vector(mode="numeric", length=0)
families_genera <- vector(mode="numeric", length=0)
for (b in 1:length(final_tip_labels_families[[a]])){ ###number_of_genera_in_family_in_tree
families_genera <- append(families_genera, unlist(strsplit(final_tip_labels_families[[a]][[b]], split="_"))[[3]])
}
families_genera <- unique(families_genera) #get the unique genera per family - will be the same as overall because there is only one per genus
families_subset_table <- table[which(table[,2]==unlist(strsplit(final_tip_labels_families[[a]][[1]], split="_"))[[2]]),] #number_of_genera_in_family_in_wcvp
wcvp_families_genera <- unique(families_subset_table[,3][[1]])
missing_genera <- length(wcvp_families_genera) - length(families_genera) ###number_of_missing_genera
message(paste(analysis_families_name[[a]], " has ", missing_genera, " missing genera", sep=""))
if (missing_genera > 0){ ###get_delay_between_families_stem_and_genera_stem_for_the_relevent_family
delay_counter <- vector(mode="numeric", length=0)
for (b in which(genera_corresponding_families_names == analysis_families_name[[a]])){ # get the stem times of the genera in the family
delay_counter <- append(delay_counter, stem_times_genera[[b]])
}
delay_mean <- mean(delay_counter)
if (length(delay_counter) > 3){
delay_sd <- sd(delay_counter)
} else {
delay_sd <- 0.1*delay_mean
}
if (nrow(families_subset_table)==1){ ###get_species_richness_of_missing_genera
species <- nrow(families_subset_table)
}
if (nrow(families_subset_table)>1){
species <- nrow(families_subset_table)-1
}
average_species_richness_of_genera_in_families <- (species-sum(number_of_species_done))/missing_genera # remaining species divided by number of missing genera
for (b in 1:missing_genera){ ###simulate_the_missing_genera
root_age <- vector(mode="numeric", length=0)
while(length(root_age) == 0){
root_age_test <- rnorm(1, delay_mean, delay_sd)
if (root_age_test > 0){
root_age <- root_age_test
}
}
missing_tree <- sim.bd.taxa.age(round(average_species_richness_of_genera_in_families)+1, 1, missing_tree_spec, missing_tree_ext, frac = 1, root_age, mrca = FALSE)[[1]] 
missing_tree <- drop.tip(missing_tree, sample(seq(1, length(missing_tree$tip.label), 1), 1))
missing_tree_table <- data.frame(missing_tree[[1]], missing_tree$edge.length, seq(1, nrow(missing_tree[[1]]), 1))
missing_tree_table <- missing_tree_table[order(missing_tree_table[,2]),]
end_times <- c(rep(0, length(missing_tree$tip.label)), branching.times(missing_tree)[-1])
start_times <- end_times + missing_tree_table[,3]
for (c in 1:length(interval_missing_diversity_storage)){
counter <- length(which(end_times <= intervals[[c]] & start_times > intervals[[c]]))
if (counter > 0){
interval_missing_diversity_storage[[c]] <- append(interval_missing_diversity_storage[[c]], counter)
}  
if ((counter == 0) & (root_age > intervals[[c]])){
interval_missing_diversity_storage[[c]] <- append(interval_missing_diversity_storage[[c]], 1)
}
}
for (c in 1:length(interval_missing_diversity_storage)){
if (length(interval_missing_diversity_storage[[c]]) > 0){
interval_missing_diversity_storage[[c]] <- sum(interval_missing_diversity_storage[[c]])
}
}
}
}
} else {
families_genera <- unlist(strsplit(tree$tip.label[[final_node_families[[a]]]], split="_"))[[3]] ###get_genera_in_family_in_tree
families_subset_table <- table[which(table[,2]==unlist(strsplit(tree$tip.label[[final_node_families[[a]]]], split="_"))[[2]]),]###get_genera_in_family_in_wcvp
wcvp_families_genera <- unique(families_subset_table[,3][[1]])
missing_genera <- length(wcvp_families_genera) - length(families_genera)###number_of_missing_genera
if (missing_genera > 0){ ###get_species_richness_of_missing_genera
if (nrow(families_subset_table)==1){
species <- nrow(families_subset_table)
}
if (nrow(families_subset_table)>1){
species <- nrow(families_subset_table)-1
}
average_species_richness_of_genera_in_families <- (species-sum(number_of_species_done))/missing_genera # remaining species divided by number of missing genera
for (b in 1:missing_genera){ ###simulate_the_missing_genera
root_age <- vector(mode="numeric", length=0)
while(length(root_age) == 0){
root_age_test <- rnorm(1, 0.5*stem_times_families[[a]], 0.25*stem_times_families[[a]]) 
if (root_age_test > 0){
root_age <- root_age_test
}
}
if (average_species_richness_of_genera_in_families >= 2){
missing_tree <- sim.bd.taxa.age(round(average_species_richness_of_genera_in_families), 1, missing_tree_spec, missing_tree_ext, frac = 1, root_age, mrca = FALSE)[[1]] 
missing_tree <- drop.tip(missing_tree, sample(seq(1, length(missing_tree$tip.label), 1), 1))
missing_tree_table <- data.frame(missing_tree[[1]], missing_tree$edge.length, seq(1, nrow(missing_tree[[1]]), 1))
missing_tree_table <- missing_tree_table[order(missing_tree_table[,2]),]
end_times <- c(rep(0, length(missing_tree$tip.label)), branching.times(missing_tree)[-1])
start_times <- end_times + missing_tree_table[,3]
for (c in 1:length(interval_missing_diversity_storage)){
counter <- length(which(end_times <= intervals[[c]] & start_times > intervals[[c]]))
if (counter > 0){
interval_missing_diversity_storage[[c]] <- append(interval_missing_diversity_storage[[c]], counter)
}  
if ((counter == 0) & (root_age > intervals[[c]])){
interval_missing_diversity_storage[[c]] <- append(interval_missing_diversity_storage[[c]], 1)
}
}
} else {
for (c in 1:length(interval_missing_diversity_storage)){
if (root_age > intervals[[c]]){
interval_missing_diversity_storage[[c]] <- append(interval_missing_diversity_storage[[c]], 1)
}
}
} 
for (c in 1:length(interval_missing_diversity_storage)){
if (length(interval_missing_diversity_storage[[c]]) > 0){
interval_missing_diversity_storage[[c]] <- sum(interval_missing_diversity_storage[[c]])
}
}
}
}
}
}
}
for (i in 1:length(interval_missing_diversity_storage)){
if (length(interval_missing_diversity_storage[[i]]) == 0){
interval_missing_diversity_storage[[i]] <- 0
}
interval_missing_diversity_rep_storage[[i]] <- append(interval_missing_diversity_rep_storage[[i]], interval_missing_diversity_storage[[i]])
}
}

####################################
###GET_TREE_DIVERSITY_AT_INTERVAL###
####################################

tree_table <- data.frame(tree[[1]], tree$edge.length)
tree_table <- tree_table[order(tree_table[,2]),] 
end_times <- c(rep(0, length(tree$tip.label)), branching.times(tree)[-1])
start_times <- end_times + tree_table[,3]
tree_table <- cbind(tree_table, end_times, start_times)

interval_diversity <- vector(mode="numeric", length=0)
for (a in 1:length(intervals)){
interval_diversity <- append(interval_diversity, length(which(end_times <= intervals[[a]] & start_times > intervals[[a]])))
}

##############################
###DIVIDE_BY_TREE_DIVERSITY###
##############################

interval_missing_diversity_rep_storage_divided <- vector("list", length(interval_missing_diversity_rep_storage))
for (a in 1:length(interval_missing_diversity_rep_storage)){
for (b in 1:length(interval_missing_diversity_rep_storage[[a]])){
if (interval_diversity[[a]] > 0){
interval_missing_diversity_rep_storage_divided[[a]] <- append(interval_missing_diversity_rep_storage_divided[[a]], (interval_diversity[[a]]/(interval_diversity[[a]]+interval_missing_diversity_rep_storage[[a]][[b]]))*100)
} else {
interval_missing_diversity_rep_storage_divided[[a]] <- append(interval_missing_diversity_rep_storage_divided[[a]], 100)
}
}
}

#####################
###GET_CI_AND_MEAN###
#####################

sd_interval_missing_lower <- vector(mode="numeric", length=0)
sd_interval_missing_upper <- vector(mode="numeric", length=0)
mean_missing <- vector(mode="numeric", length=0)

for (i in 1:length(interval_missing_diversity_rep_storage_divided)){
mean_missing <- append(mean_missing, mean(interval_missing_diversity_rep_storage_divided[[i]]))
sd_interval_missing_lower <- append(sd_interval_missing_lower, mean_missing[[i]] - sd(interval_missing_diversity_rep_storage_divided[[i]]))
sd_interval_missing_upper <- append(sd_interval_missing_upper, mean_missing[[i]] + sd(interval_missing_diversity_rep_storage_divided[[i]]))
}

for (i in 1:length(interval_missing_diversity_rep_storage_divided)){
if (is.na(mean_missing[[i]]) == TRUE){
mean_missing[[i]] <- 0
}
if (is.na(sd_interval_missing_lower[[i]]) == TRUE){
sd_interval_missing_lower[[i]] <- 0
}
if (is.na(sd_interval_missing_upper[[i]]) == TRUE){
sd_interval_missing_upper[[i]] <- 0
}
}

############################
###COMPILATION_INTO_TABLE###
############################

final_table <- data.frame(intervals, mean_missing, sd_interval_missing_lower, sd_interval_missing_upper)
write.table(final_table, "0.99_sampling_through_time.tsv", row.names=FALSE, col.names=c("Time", "Mean", "SD_low", "SD_up"), quote=FALSE)
##########
###PLOT###
##########

final_plot <- ggplot(data=final_table, aes(x = final_table[,1], y = final_table[,2]))+
geom_line(aes(y=final_table[,2]), colour="black", size=0.5) + 
geom_ribbon(aes(ymin=final_table[,3], ymax=final_table[,4]),size=0.2,alpha=0.15,fill="black")+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(size=1,colour="black"), text=element_text(size=16,colour="black"), axis.ticks=element_line(size=1,colour="black"))+
labs(y="% incorporated", x="Ma")+
scale_x_reverse()
