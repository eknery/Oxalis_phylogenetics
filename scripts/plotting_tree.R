
BiocManager::install("ggtree")
if(!require("ggtree")) install.packages("ggtree"); library("ggtree")
if(!require("treeio")) install.packages("treeio"); library("treeio")
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")

### input directory
dir_input = "8_bayesian_phylogenies/"

### read .nexus tree
beast_tree <- read.beast(
  paste0(dir_input, "mcc_median_annotated.tre")
)

############################ PROCESSING TREE ##################################

### posterior limits
beast_tree@data$posterior_lim = beast_tree@data$posterior
beast_tree@data$posterior_lim[beast_tree@data$posterior >= 0.95] = "*"
beast_tree@data$posterior_lim[beast_tree@data$posterior < 0.95] = NA

### removing specimens names tip labels
## get current names
tip_label = beast_tree@phylo$tip.label
## split
splited = strsplit(tip_label, split = "_")
## loop over splited names
new_label = c()
for(i in 1:length(splited)){
  last_index = length(splited[[i]])
  newname = paste0(splited[[i]][-last_index], collapse = "_")
  newname = gsub(pattern= "'", replacement = "", x = newname)
  new_label = c(new_label, newname)
}
## remove "_"
new_label = gsub("_"," ",new_label)
## replace tip labels
beast_tree@phylo$tip.label = new_label

###################################### PLOT TREE ###############################

options(ignore.negative.edge=TRUE)

### export pie chart
tiff("figures/phylo_plot.tiff",
     width= 28.5, 
     height= 26,
     units="cm",
     res=900)
ggtree(
  beast_tree,
  ladderize = T,
  right = T,
  linewidth = 0.3,
  aes(x = length),
  ) +
  geom_range(
    range = 'height_0.95_HPD', 
    color = 'gray', 
    alpha = 0.5, 
    size = 0.7
  ) +
  geom_nodelab(
    aes(x=branch, label=posterior_lim),
    vjust=0.3, 
    size=1.5
  ) +
  geom_tiplab(
    size = 0.7,
    offset = 0,
    align = F
  ) +
  scale_x_continuous(
    breaks =  seq(0, 100, by = 10),
    labels =  as.character(seq(100, 0, by = -10))
  ) +
  theme_tree2(
    bgcolor = "white",
    fgcolor = "black"
  ) 
dev.off()

