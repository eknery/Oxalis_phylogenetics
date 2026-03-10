### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### input diretory
dir_input = "TRNL-F/"
file_names = list.files(dir_input)

### loading data
all_seqs = c()
for(i in 1:length(file_names) ){
  seq = readLines(paste0(dir_input, file_names[i]))
  seq_name = paste0(">", str_remove(string = file_names[i], pattern = ".seq"))
  one_seq = paste(seq_name, seq, sep = "\n")
  all_seqs = paste0(all_seqs, "\n", one_seq)
}

writeLines(all_seqs,"all_seqs.txt")


