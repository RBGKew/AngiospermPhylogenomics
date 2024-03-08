clades_colours_long_label <- list(Gymnosperms = "black",
                       Angiosperms = "#882255",
                       Ceratophyllales = "#888888",
                       Magnoliids = "#117733",
                       Monocots = "#DDCC77",
                       Monocots_Commelinids = "#999933",
                       Eudicots = "#332288",
                       Eudicots_Rosids = "#88CCEE",
                       Eudicots_Rosids_Malvidae = "#44AA99",
                       Eudicots_Rosids_Fabidae = "#6699CC",
                       Eudicots_Asterids = "#CC6677",
                       Eudicots_Asterids_Campanuliidae = "#AA4499",
                       Eudicots_Asterids_Lamiidae = "#661100")

clades_colours <- clades_colours_long_label
names(clades_colours) <- unlist(lapply(names(clades_colours_long_label),
							FUN = function(x) tail(unlist(strsplit(x, "_")), n = 1)))
