### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose directory with sequences
dir_input = "1_selected_accessions/" 

### choose directory with sequences
dir_out = "1_selected_accessions/"

### read acessions table
best_acc = read.csv(paste0(dir_input, "best_acc.csv") )

### only Oxalis 
best_acc_ox = best_acc[grepl("Oxalis_", best_acc$taxon),]

sum(!is.na(best_acc_ox$ITS))
sum(!is.na(best_acc_ox$ncpGS))
sum(!is.na(best_acc_ox$psbA_trnH))
sum(!is.na(best_acc_ox$psbJ_petA))
sum(!is.na(best_acc_ox$rbcL))
sum(!is.na(best_acc_ox$trnL_trnF))
sum(!is.na(best_acc_ox$trnS_trnG))
sum(!is.na(best_acc_ox$trnT_trnL))
