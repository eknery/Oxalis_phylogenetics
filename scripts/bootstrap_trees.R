### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### input diretory
dir_input = "4_trimmed_sequences/"
### available files
list.files(dir_input)
### choose a file
file_name = "ITS.fasta"

### loading data
fasta_file = read.phyDat(
  paste0(dir_input, file_name),
  format = "fasta",
  type = "DNA"
  )

### locus name
locus_name = str_remove(
  string = file_name, 
  pattern = ".fasta"
  )

############################## INFERING ML TREE ###############################

### ML fits
ml_fits = bootstrap.phyDat(
  x = fasta_file, 
  FUN = function(x)optim.pml(pml(NJ(dist.ml(x)), data = x), model = "GTR"), 
  bs = 100
)

### ML trees
ml_trees = list()
for(i in 1:length(ml_fits)){
  ml_trees[[i]] = ml_fits[[i]]$tree
}

### export trees
write.tree(
  phy =  ml_trees,
  file = paste0("7_ml_phylogenies/bootstrap/", locus_name,".tre")
  )


