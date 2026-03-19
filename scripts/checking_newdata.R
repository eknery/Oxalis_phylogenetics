### libraries 
if(!require("ape")) install.packages("ape"); library("ape")

### input directory
dir_input = "1_selected_accessions/"

### read acessions table
best_acc = read.csv(paste0(dir_input, "best_acc.csv") )

############################## PROCESSING DATA #################################

### loci names
all_loci = colnames(best_acc)[!colnames(best_acc) %in% c("taxon", 
                                                         "specimen",
                                                         "Nmarker",
                                                         "Kmarker",
                                                         "duplicate")]
### choose a locus name
locus_name = all_loci[9]
### lines with information
locus_acc = best_acc[!is.na(best_acc[,locus_name]),]
### unpublished accessions
unp_acc = locus_acc[locus_acc[,locus_name] == 1,]
unp_names = paste0(unp_acc[,"taxon"], "_", unp_acc[,"specimen"])
### read newdata
newdata = read.FASTA(
  paste0("0_newdata/",locus_name,".fasta")
)
### selected accessiosn not in newdata
not_in_data = unp_names[!unp_names %in% names(newdata)]
###
not_in_data
