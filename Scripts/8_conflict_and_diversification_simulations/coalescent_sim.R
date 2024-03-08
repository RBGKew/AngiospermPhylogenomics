library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(TreePar)
library(ggplot2)

###params###
#intermediate burst#
#speciation only lambda and 250 tips = 0.02, 0.2, 0.02, mu = 0.0, 0.0, 0.0#
#birth death 
#extinction drop lambda and 250 tips = 0.2, 0.2, 0.2, mu = 0.18, 0, 0.18

####################
###INITIAL_VALUES###
####################

generation_time <- 50
n_gene_trees <- 100
n_reps <- 50
effective_population_size_approx <- 0.005 # adjust this value as set out in methods. to translate to a scaled value, multiply by 1000000
n_tips <- 100
lambda <- c(0.75, 0.75, 0.75, 0.75) # adjust these rates, as set out in supplementary results figure and methods
mu <- c(0, 0.675, 0.675, 0)
frac <- c(1, 1, 1, 1)
times <- c(0, 2, 4, 6)
K <- 0

window = 0.5
conflict_slices <- seq(0, 200, window)
slice_conflict_values <- vector("list", length(conflict_slices))
slice_plot_times <- vector(mode="numeric", length=0)
interval_conflict_values <- vector("list", length(conflict_slices)-1)
interval_plot_times <- vector(mode="numeric", length=0)

####################
###RUN_SIMULATION###
####################

for(k in 1:n_reps){

###SIMULATE_SPECIES_TREE###
###########################

test_age <- 0
while (test_age <= 8){
entire_tree <- sim.rateshift.taxa(n_tips, n_reps, lambda, mu, frac, times, complete = TRUE, K=0)
entire_tree <- entire_tree[[1]]
entire_tree <- extract.clade(entire_tree, findMRCA(entire_tree, drop.tip(entire_tree, getExtinct(entire_tree))$tip.label, "node"))
test_age <- max(node.depth.edgelength(drop.tip(entire_tree, getExtinct(entire_tree))))
}

sampling_times <- round(node.depth.edgelength(entire_tree)[seq(1, length(entire_tree$tip.label), 1)] - max(node.depth.edgelength(entire_tree)), digits = 5)	# tip sampling times for incorporation into gene trees

###DEFINE_OTHER_USEFUL_VARIABLES###
###################################

tip_vector <- seq(1, length(entire_tree$tip.label), 1)
node_numbers <- seq(length(tip_vector)+1, length(tip_vector) + (length(tip_vector)-1), 1) 
file_names <- paste(seq(1, n_gene_trees, 1), ".tre", sep="")

###DEFINE_COALESCENT_POPULATION_CONSTRAINTS###
##############################################

###define the actual populations - things that can concievably join###

coalescent_populations <- vector("list", (length(tip_vector)-1))
for (i in 1:length(coalescent_populations)){
coalescent_populations[[i]] <- extract.clade(entire_tree, node_numbers[[i]])$tip.label
}
coalescent_populations <- coalescent_populations[order(sapply(coalescent_populations,length),decreasing=F)]

###define the minimum and maximum ages for the populations###

coalescent_populations_minimum_ages <- vector("list", length(coalescent_populations))
coalescent_populations_maximum_ages <- vector("list", length(coalescent_populations))
for (i in 1:length(coalescent_populations_minimum_ages)){
coalescent_populations_minimum_ages[[i]] <- max(node.depth.edgelength(entire_tree)) - findMRCA(entire_tree, coalescent_populations[[i]], "height")	#age of the population in the species tree
if(length(coalescent_populations[[i]]) < length(entire_tree$tip.label)){
coalescent_populations_maximum_ages[[i]] <- max(node.depth.edgelength(entire_tree)) - findMRCA(entire_tree, extract.clade(entire_tree, entire_tree[[1]][,1][[which(entire_tree[[1]][,2] == findMRCA(entire_tree, coalescent_populations[[i]], "node"))]])$tip.label, "height")	#age of ancestral population in species tree
}																																																																				#NB no maximum age of the coalescent population at the root
}

###SIMULATE_GENE_TREES###
#########################

###start a list of gene trees###

gene_trees <- vector("list", n_gene_trees)
for (y in 1:length(gene_trees)){ 
message(paste("gene tree ", y, sep=""))

###generate the starting pool for a given gene tree###

pool <- vector("list", length(entire_tree$tip.label))
for (i in 1:length(pool)){
pool[[i]] <- list(edge=matrix(c(2,1),1,2), tip.label=entire_tree$tip.label[[i]], edge.length=0.0, Nnode=1)
class(pool[[i]]) <- "phylo"
}
message("initial pool defined")

###start going through the coalescent populations###

checkpoint <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)
for (i in 1:length(coalescent_populations_minimum_ages)){											# large step starts here through each coalescent population
if (i %in% checkpoint){
message(paste(i, " coalescent populations done", sep=""))
}

a=0
continue_at_population <- "yes"
while(continue_at_population == "yes"){
a=a+1
if (length(pool) >= 2){															# if theres still things to join

###subsection to target the search and speed up the simulation###

pools_to_sample <- vector(mode="numeric", length=0)
for (b in 1:length(pool)){  
counter <- length(which(pool[[b]]$tip.label %in% coalescent_populations[[i]]))	# are all the pools tips in the coalescent population
if (counter == length(pool[[b]]$tip.label)){
pools_to_sample <- append(pools_to_sample, b)
}
}

if (length(pools_to_sample) > 1){

combinations <- vector("list", 0)
for (b in 1:length(pools_to_sample)){
for (c in seq(1, length(pools_to_sample), 1)[-b]){
combinations[[length(combinations)+1]] <- sort(c(pools_to_sample[[b]], pools_to_sample[[c]]))
}
}
combinations <- unique(combinations)

combinations_ages <- vector(mode="numeric", length=0)
for (b in 1:length(combinations)){
if ((length(pool[[combinations[[b]][[1]]]]$tip.label) == 1) & (length(pool[[combinations[[b]][[2]]]]$tip.label) == 1)){
combinations_ages <- append(combinations_ages, (rexp(1, 1/(effective_population_size_approx))*2*generation_time) + coalescent_populations_minimum_ages[[i]])
}
if ((length(pool[[combinations[[b]][[1]]]]$tip.label) > 1) & (length(pool[[combinations[[b]][[2]]]]$tip.label) == 1)){
if ((node.depth.edgelength(pool[[combinations[[b]][[1]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[1]]]]$tip.label[[1]])]])) < coalescent_populations_minimum_ages[[i]]){
combinations_ages <- append(combinations_ages, (rexp(1, 1/(effective_population_size_approx))*2*generation_time) + coalescent_populations_minimum_ages[[i]])
} else {
combinations_ages <- append(combinations_ages, (rexp(1, 1/(effective_population_size_approx))*2*generation_time) + (node.depth.edgelength(pool[[combinations[[b]][[1]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[1]]]]$tip.label[[1]])]])))
}
}
if ((length(pool[[combinations[[b]][[2]]]]$tip.label) > 1) & (length(pool[[combinations[[b]][[1]]]]$tip.label) == 1)){
if ((node.depth.edgelength(pool[[combinations[[b]][[2]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[2]]]]$tip.label[[1]])]])) < coalescent_populations_minimum_ages[[i]]){
combinations_ages <- append(combinations_ages, (rexp(1, 1/(effective_population_size_approx))*2*generation_time) + coalescent_populations_minimum_ages[[i]])
} else {
combinations_ages <- append(combinations_ages,(rexp(1, 1/(effective_population_size_approx))*2*generation_time) + (node.depth.edgelength(pool[[combinations[[b]][[2]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[2]]]]$tip.label[[1]])]])))
}
}
if ((length(pool[[combinations[[b]][[2]]]]$tip.label) > 1) & (length(pool[[combinations[[b]][[1]]]]$tip.label) > 1)){
if (((node.depth.edgelength(pool[[combinations[[b]][[1]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[1]]]]$tip.label[[1]])]])) < coalescent_populations_minimum_ages[[i]]) & ((node.depth.edgelength(pool[[combinations[[b]][[2]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[2]]]]$tip.label[[1]])]])) < coalescent_populations_minimum_ages[[i]])){
combinations_ages <- append(combinations_ages, (rexp(1, 1/(effective_population_size_approx))*2*generation_time) + coalescent_populations_minimum_ages[[i]])
} else if ((node.depth.edgelength(pool[[combinations[[b]][[1]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[1]]]]$tip.label[[1]])]])) > (node.depth.edgelength(pool[[combinations[[b]][[2]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[2]]]]$tip.label[[1]])]]))){ 
combinations_ages <- append(combinations_ages,(rexp(1, 1/(effective_population_size_approx))*2*generation_time) + (node.depth.edgelength(pool[[combinations[[b]][[1]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[1]]]]$tip.label[[1]])]])))
} else {
combinations_ages <- append(combinations_ages,(rexp(1, 1/(effective_population_size_approx))*2*generation_time) + (node.depth.edgelength(pool[[combinations[[b]][[2]]]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == pool[[combinations[[b]][[2]]]]$tip.label[[1]])]])))
}
}
}

go <- "no"
if (i < length(coalescent_populations_maximum_ages)){
if (min(combinations_ages) < coalescent_populations_maximum_ages[[i]]){
go <- "go"
}
} else {
go <- "go"
}

if (go == "go"){

to_join <- pool[combinations[[which(combinations_ages == min(combinations_ages))]]]
to_join_age <- min(combinations_ages)

###join_branches###

if ((length(to_join[[1]]$tip.label) == 1) & (length(to_join[[2]]$tip.label) == 1)){
to_join[[1]]$edge.length <- to_join_age - abs(sampling_times)[[which(entire_tree$tip.label == to_join[[1]]$tip.label)]]
to_join[[2]]$edge.length <- to_join_age - abs(sampling_times)[[which(entire_tree$tip.label == to_join[[2]]$tip.label)]]
}
if ((length(to_join[[1]]$tip.label) > 1) & (length(to_join[[2]]$tip.label) == 1)){
to_join[[1]] <- addroot(to_join[[1]], to_join_age - (node.depth.edgelength(to_join[[1]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == to_join[[1]]$tip.label[[1]])]])))
to_join[[2]]$edge.length <- to_join_age - abs(sampling_times[[which(entire_tree$tip.label == to_join[[2]]$tip.label)]])
}
if ((length(to_join[[1]]$tip.label) == 1) & (length(to_join[[2]]$tip.label) > 1)){
to_join[[1]]$edge.length <- to_join_age - abs(sampling_times[[which(entire_tree$tip.label == to_join[[1]]$tip.label)]])
to_join[[2]] <- addroot(to_join[[2]], to_join_age - (node.depth.edgelength(to_join[[2]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == to_join[[2]]$tip.label[[1]])]])))
}
if ((length(to_join[[1]]$tip.label) > 1) & (length(to_join[[2]]$tip.label) > 1)){
to_join[[1]] <- addroot(to_join[[1]], to_join_age - (node.depth.edgelength(to_join[[1]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == to_join[[1]]$tip.label[[1]])]])))
to_join[[2]] <- addroot(to_join[[2]], to_join_age - (node.depth.edgelength(to_join[[2]])[[1]] - abs(sampling_times[[which(entire_tree$tip.label == to_join[[2]]$tip.label[[1]])]])))
}

new <- bind.tree(to_join[[1]], to_join[[2]])

pool <- pool[-combinations[[which(combinations_ages == to_join_age)]]]    # the pool is updated internally therefore will search for the pool elements to join on the new pool
pool[[length(pool)+1]] <- new
} else {
continue_at_population <- "no"
}
} else {
continue_at_population <- "no"
}
} else {
continue_at_population <- "no"
}
}
}
gene_trees[[y]] <- pool[[1]]
message("gene tree done")
}

##################################
###COMPILE_CONFLICT_INFORMATION###
##################################

###get extant trees###
######################

entire_tree_extant <- drop.tip(entire_tree, getExtinct(entire_tree))

gene_trees_extant <- vector("list", length(gene_trees))
for (i in 1:length(gene_trees)){
gene_trees_extant[[i]] <- drop.tip(gene_trees[[i]], getExtinct(entire_tree))
}

species_tree_clades <- vector("list", length(entire_tree_extant$tip.label) - 2)
for (i in 1:length(species_tree_clades)){
species_tree_clades[[i]] <- extract.clade(entire_tree_extant, length(entire_tree_extant$tip.label)+1+i)
}

###get amount of conflict on different clades###
################################################

tip_vector <- seq(1, length(entire_tree_extant$tip.label), 1)	# get vector of tips in species tree

###get species tree clades###

strict_support <- vector("list", 0)
node <- vector("list", 0)
for (a in seq((length(tip_vector)+2), (length(tip_vector) + (length(tip_vector)-1)), 1)){				# for every descendant node in the tree
strict_support_vector <- vector(mode="numeric", length=0)												# number of gene trees that support the branch										# number of gene trees that support the branch
tested_vector <- vector(mode="numeric", length=0)

message(paste("testing branch ", a, sep=""))
overall_descendant_clade <- extract.clade(entire_tree_extant, a)														# the descendant clade of the branch of interest

overall_descendant_clade_tips <- overall_descendant_clade$tip.label
overall_other_tips <- entire_tree_extant$tip.label[-which(entire_tree_extant$tip.label %in% overall_descendant_clade_tips)] 

###now run the check against the gene trees###

for (b in 1:length(gene_trees_extant)){
if (length(which(overall_descendant_clade_tips %in% gene_trees_extant[[b]]$tip.label)) > 0){
if (length(which(overall_other_tips %in% gene_trees_extant[[b]]$tip.label)) > 0){
tested_vector <- append(tested_vector, b)
if (length(overall_descendant_clade_tips[which(overall_descendant_clade_tips %in% gene_trees_extant[[b]]$tip.label)]) > 1){
gene_tree_descendant_clade <- extract.clade(gene_trees_extant[[b]], findMRCA(gene_trees_extant[[b]], overall_descendant_clade_tips[which(overall_descendant_clade_tips %in% gene_trees_extant[[b]]$tip.label)], "node"))$tip.label
} else {
gene_tree_descendant_clade <- overall_descendant_clade_tips[which(overall_descendant_clade_tips %in% gene_trees_extant[[b]]$tip.label)]
}
gene_tree_other_tips <- gene_trees_extant[[b]]$tip.label[-which(gene_trees_extant[[b]]$tip.label %in% gene_tree_descendant_clade)]
if (length(which(gene_tree_descendant_clade %in% overall_descendant_clade_tips[which(overall_descendant_clade_tips %in% gene_trees_extant[[b]]$tip.label)])) == length(gene_tree_descendant_clade)){
if (length(which(gene_tree_other_tips %in% overall_other_tips[which(overall_other_tips %in% gene_trees_extant[[b]]$tip.label)])) == length(gene_tree_other_tips)){
strict_support_vector <- append(strict_support_vector, b)
}
}
}
}
}
strict_support[[length(strict_support) + 1]] <- length(strict_support_vector)/length(tested_vector)
node[[length(node) + 1]] <- a
}

strict_support <- append(unlist(strict_support), rep(1, length(entire_tree_extant$tip.label)), after=0)

################################
###compile in to a data frame###
################################

conflict_data_frame <- data.frame(entire_tree_extant[[1]], entire_tree_extant$edge.length)													#edge length - i=3
conflict_data_frame <- conflict_data_frame[order(conflict_data_frame[,2]),]
conflict_data_frame <- cbind(conflict_data_frame, strict_support)														#clade support - i=4
conflict_data_frame <- cbind(conflict_data_frame, c(rep(0, length(entire_tree_extant$tip.label)), branching.times(entire_tree_extant)[-1]))	#add branch end times - i=5
conflict_data_frame <- cbind(conflict_data_frame, as.numeric(conflict_data_frame[,3]) + as.numeric(conflict_data_frame[,5]))				#add branch start times - i=6

#####################################
###compile list by slice or window###
#####################################

for (a in 1:length(slice_conflict_values)){
slice_conflict_values[[a]] <- append(slice_conflict_values[[a]], conflict_data_frame[,4][which(conflict_data_frame[,5] <= conflict_slices[[a]] & conflict_data_frame[,6] > conflict_slices[[a]])])
if (a < length(slice_conflict_values)){
slice_plot_times <- append(slice_plot_times, c(conflict_slices[[a]], conflict_slices[[a+1]]-0.0001))
} else {
slice_plot_times <- append(slice_plot_times, c(conflict_slices[[a]], conflict_slices[[a]]+0.001))
}
}

for (a in 2:length(conflict_slices)){
interval_conflict_values[[a-1]] <- append(interval_conflict_values[[a-1]], conflict_data_frame[,4][which(conflict_data_frame[,5] <= conflict_slices[[a]] & conflict_data_frame[,6] >= conflict_slices[[a-1]])])
interval_plot_times <- append(interval_plot_times, c(conflict_slices[[a-1]], conflict_slices[[a]]-0.0001))
}

###################
###WRITE_TO_FILE###
###################

write.table(conflict_data_frame, paste("med_pop/", k, "conflict.tsv", sep=""), col.names=FALSE, quote=FALSE)
write.tree(entire_tree, paste("med_pop/", k, "species_tree.tre", sep=""))
write.tree(entire_tree_extant, paste("med_pop/", k,"species_tree_extant.tre", sep=""))

for (i in 1:length(file_names)){
write.tree(gene_trees[[i]], paste("med_pop/", k, "_", file_names[[i]], sep=""))
}

for (i in 1:length(file_names)){
write.tree(gene_trees_extant[[i]], paste("med_pop/", k, "_extant_", file_names[[i]], sep=""))
}

}

###prepare plot###
##################

#interval_conflict_values_plot <- unlist(lapply(interval_conflict_values, mean))

slice_conflict_values_plot <- unlist(lapply(slice_conflict_values, mean))

#conflict_data_frame_final <- data.frame(interval_plot_times, rep(interval_conflict_values_plot, each=2))  
conflict_data_frame_final <- data.frame(slice_plot_times, (1-rep(slice_conflict_values_plot, each=2))*100)  

conflict_plot <- ggplot(data=conflict_data_frame_final, aes(x=conflict_data_frame_final[,1], y=conflict_data_frame_final[,2])) +
geom_line(aes(y=conflict_data_frame_final[,2]))+
geom_rect(xmin=9, ymin=0, xmax=11, ymax=50, fill="black", alpha=0.5)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(size=1,colour="black"), text=element_text(size=16,colour="black"), axis.ticks=element_line(size=1,colour="black"))+
expand_limits(x = 0, y = 0) +
scale_y_continuous(expand = c(0, 0), limits=c(0, 20)) +
scale_x_reverse(expand=c(0,0), limits = c(30, 0))+
labs(y="% incongruence", x="Time(Myr)")



