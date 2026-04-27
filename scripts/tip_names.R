### libraries
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ips")) install.packages("ips"); library("ips")

### file names
dir_input = "8_bayesian_phylogenies/"
tree_name = "mcc.tre"

### import
tree = read.tree(file = paste0(dir_input, tree_name))

### check tree
plot(tree, cex=0.1)

### tip labels
tip_label = tree$tip.label

### split
splited = strsplit(tip_label, split = "_")

### loop over splited names
new_label = c()
for(i in 1:length(splited)){
  last_index = length(splited[[i]])
  newname = paste0(splited[[i]][-last_index], collapse = "_")
  newname = gsub(pattern= "'", replacement = "", x = newname)
  new_label = c(new_label, newname)
}

### replace tip labels
newtree = tree
newtree$tip.label = new_label

### check tree
plot(newtree, cex=0.1)

### export 
write.tree(
  phy =  newtree,
  file = paste0(dir_input,"mcc_clean.tre")
)
