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

### sequence per locus
sum(!is.na(acc_ox$ITS))
sum(!is.na(acc_ox$ncpGS))
sum(!is.na(acc_ox$psbA_trnH))
sum(!is.na(acc_ox$psbJ_petA))
sum(!is.na(acc_ox$rbcL))
sum(!is.na(acc_ox$trnL_trnF))
sum(!is.na(acc_ox$trnS_trnG))
sum(!is.na(acc_ox$trnT_trnL))
