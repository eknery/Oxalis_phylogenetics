### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### input diretory
dir_input = "4_trimmed_sequences/"
file_names = list.files(dir_input)

### loading data
fasta_list = list()
for(i in 1:length(file_names) ){
  fasta_list[[i]] = read.phyDat(paste0(dir_input, file_names[i]),
                                format = "fasta",
                                type = "DNA"
  )
  names(fasta_list)[i] = str_remove(string = file_names[i], 
                                    pattern = ".fasta")
}

############################## INFERING ML TREE ###############################

## select output dir
dir_out = "7_ml_phylogenies/bootstrap/"

### loci names
loci_names = sub(".fasta", "", file_names)

### loop
for(locus_name in loci_names){
  ### pick one aligment
  one_fasta = fasta_list[[locus_name]]
  ### ML fits
  ml_fits = bootstrap.phyDat(
    x = one_fasta, 
    FUN = function(x)optim.pml(pml(NJ(dist.ml(x)), data = x), model = "GTR"), 
    bs = 100
  )
  ### ML trees
  ml_trees = list()
  for(i in 1:length(ml_fits)){
    ml_trees[[i]] = ml_fits[[i]]$tree
  }
  ## export trees
  write.tree(
    phy =  ml_trees,
    file = paste0(dir_out,locus_name,".tre")
  )
  print(paste0("ML bootstrap done: ", locus_name))
}

