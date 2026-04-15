### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### input diretory
dir_input = "4_trimmed_sequences/"
file_name = "ITS.fasta"

### read fasta
one_fasta = read.phyDat(
  paste0(dir_input, file_name),
  format = "fasta",
  type = "DNA"
)  
## locus name
locus_name = sub(".fasta", "", file_name)
## number of sequences
n = length(one_fasta)
## number of sequences from focal group
nfocal = sum(grepl("Oxalis_", names(one_fasta) ))
## sequence length
seqLength = length(one_fasta[[1]])
## fit models
modelfits = modelTest(
  object = one_fasta, 
  model = "all", 
  G = TRUE, 
  k = 4, 
  I = TRUE
)
## best fit model
bestfit = modelfits[modelfits$AICcw == max(modelfits$AICcw),]

#### data frame for results
best_models = data.frame(
  locus_name = locus_name, 
  n = n,
  nfocal = nfocal,
  seqLength = seqLength,
  bestModel = bestfit$Model,
  stringsAsFactors = FALSE
)

### export dataframe
write.csv(
  best_models, 
  paste0("6_subs_models/",locus_name, ".csv"),
  row.names = FALSE
)
