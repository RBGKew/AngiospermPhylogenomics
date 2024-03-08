library(BAMMtools)
library(phytools)

tree <- read.tree("ladderised.tre")
#tree <- ladderize(tree, right=TRUE)
event_data <- getEventData(tree, eventdata="div_all.txt", burnin=0.5, nsamples=2000, type="diversification")

best <- getBestShiftConfiguration(event_data, expectedNumberOfShifts=10, threshold=4)

###get_nodes###

table <- read.table("New_Calibrations.csv")

###get tips and node definitions down to only genus level bc we are working with a one per genus tree, calibrations were defined within the original pre pruned calibration tree###

genus_only_tree_tips <- unlist(lapply(lapply(lapply(lapply(tree$tip.label, strsplit, "_"), unlist), "[", c(1,2,3)), paste, collapse="_"))
genus_only_calibration_tips_one <- unlist(lapply(lapply(lapply(lapply(table[,1], strsplit, "_"), unlist), "[", c(1,2,3)), paste, collapse="_"))
genus_only_calibration_tips_two <- unlist(lapply(lapply(lapply(lapply(table[,2], strsplit, "_"), unlist), "[", c(1,2,3)), paste, collapse="_"))

nodes <- vector(mode="numeric", length=0)
for (i in 1:nrow(table)){
nodes <- append(nodes, findMRCA(tree, c(tree$tip.label[[which(genus_only_tree_tips == genus_only_calibration_tips_one[[i]])]], tree$tip.label[[which(genus_only_tree_tips == genus_only_calibration_tips_two[[i]])]]), "node"))
}

###adjust colour.interval for the colour interval issue we have been discussing 

pdf("circle_bamm_prior10_net_div_ladder.pdf", height=200, width=200)
plot(best, method="polar", labels = TRUE, cex=0.225, lwd=0.8, legend=TRUE, spex="netdiv", color.interval=c(0, 0.4))
nodelabels(node = nodes, pch=18, col="black", cex=5)
dev.off()





