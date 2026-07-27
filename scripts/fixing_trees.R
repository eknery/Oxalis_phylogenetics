### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### directory
dir_input = "7_ml_phylogenies/8_loci_434_spp/"

### files
bestml = read.tree(file = paste0(dir_input, "bestMl.tre"))

boots = read.tree(file = paste0(dir_input, "bootstrap.tre"))
  
boots

### tips to drop
to_drop = c(
  "Biophytum_abyssinicum_BIOAB",
  "Platytheca_galioides_Crayn701c",
  "Tetratheca_hirsuta_Butcher915",
  "Tremandra_diffusa_Butcher961"
)

###
bestml_pruned = drop.tip(
  phy = bestml, 
  tip = to_drop
  )

### prunning
boots_pruned = boots
for(i in 1:length(boots)){
  boots_pruned[[i]] = drop.tip(
    phy = boots[[i]], 
    tip = to_drop
  )
}

bestml_supp = plotBS(
  tree = bestml_pruned, 
  trees =  boots_pruned, 
  type = "phylogram", 
  method = "FBP"
)

bestml_supp$node.label[is.na(bestml_supp$node.label)] = 0 
bestml_supp$node.label = bestml_supp$node.label * 100

dir_out = "7_ml_phylogenies/8_loci_430_spp/"
write.tree(bestml_pruned, file = paste0(dir_out, "bestml.tre") )
write.tree(boots_pruned, file = paste0(dir_out, "bootstrap.tre") )
write.tree(bestml_supp, file = paste0(dir_out, "bestML_support.tre") )

