### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose directory with sequences
dir_input = "1_selected_accessions/" 

### read acessions table
acc = read.csv(paste0(dir_input, "best_acc.csv") )

### only Oxalis 
acc_ox = acc[grepl("Oxalis_", acc$taxon),]

### species level
sum(!grepl("aff_|subsp_|var_|forma_", unique(acc_ox$taxon) ))
### other taxa
sum(grepl("aff_|subsp_|var_|forma_", unique(acc_ox$taxon) ))

### data per column
x = apply(acc_ox, FUN = function(x) {sum(!is.na(x))} , MARGIN = 2)

### total sequences
sum(x[c(-1,-2)])

