### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ips")) install.packages("ips"); library("ips")

### file names
dir_input = "7_ml_phylogenies/raxml - complete_8_loci_421_spp/"
tree_name = "bestTree_support_unrooted.tre"

### import
tree = read.tree(file = paste0(dir_input, tree_name))
plot(tree, cex=0.1)

### cutoff for collaping
cutoff = 75
### colapse!
ctree = collapseUnsupportedEdges(
  phy = tree,
  value = "node.label", 
  cutoff = cutoff
)

### plot
plot(
  ctree, 
  cex =0.1
)

### get rid of info
ctree$edge.length = NULL
ctree$node.label = NULL

### export
write.tree(
  phy =  ctree,
  file = paste0(dir_input,"bestTree_topology", cutoff,"_unrooted.tre")
)
           
