### library
if(!require("phylotools")) install.packages("phylotools"); library("phylotools")


### chose input directory
dir_input = "5_concatenated_sequences/"
### choose output directory
dir_out = "5_concatenated_sequences/"

### read FASTA
fasta <- read.fasta(
  paste0(dir_input, "nuclear_2_loci_396_spp.fasta")
)

### export
write.dna(
  x = fasta, 
  file = paste(dir_out, "align.phy"), 
  format = "interleaved"
)
