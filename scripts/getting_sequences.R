### libraries 
if(!require("ape")) install.packages("ape"); library("ape")

### input directory
dir_input = "1_selected_accessions/"

### exrpoting directory
dir_out = "2_raw_sequences/"

### read acessions table
best_acc = read.csv(paste0(dir_input, "best_acc2.csv") )

############################### GETTING SEQUENCES #############################

### loci names
all_loci = colnames(best_acc)[!colnames(best_acc) %in% c("taxon", 
                                                         "specimen",
                                                         "Nmarker",
                                                         "Kmarker",
                                                         "duplicate")]
### newdata files
file_names = list.files("0_newdata")

### looping over loci
for(locus_name in all_loci){
  ### lines with information
  locus_acc = best_acc[!is.na(best_acc[,locus_name]),]
  ### published accessions, numbers, and names
  pub_acc = locus_acc[locus_acc[,locus_name] != 1,]
  pub_nums = pub_acc[,locus_name]
  pub_names = paste0(pub_acc[,"taxon"], "_", pub_acc[,"specimen"])
  if(length(pub_nums) > 0) {
    ### downloading published sequences
    pub_seqs = read.GenBank( access.nb = pub_nums )
    ### naming sequences
    names(pub_seqs) = pub_names
    ### unpublished accessions
    unp_acc = locus_acc[locus_acc[,locus_name] == 1,]
    unp_names = paste0(unp_acc[,"taxon"], "_", unp_acc[,"specimen"])
  } else {
    pub_seqs = c()
  }
  if( paste0(locus_name,".fasta") %in% file_names){
    ### read newdata
    newdata = read.FASTA(
      paste0("0_newdata/",locus_name,".fasta")
    )
  } else {
    newdata = c()
  }
  if(length(newdata) > 0){
    ### select only selected accessions in newdata
    unp_seqs = newdata[names(newdata) %in% unp_names]
  } else {
    unp_seqs = c()
  }
  ### joining published and unpublished sequences
  locus_seqs = c(pub_seqs, unp_seqs)
  ### ordering
  locus_seqs = locus_seqs[sort(names(locus_seqs))]
  ### export
  write.dna(
    locus_seqs, 
    file = paste0(dir_out, locus_name, ".fasta"), 
    format = 'fasta' 
  )
  ### check
  print(paste0("Locus complete: ", locus_name))
}

