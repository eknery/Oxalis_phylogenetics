### library
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### choose directory with sequences
dir_input = "5_concatenated_sequences/"

### list file names
all_files = list.files(path = paste0(dir_input), pattern = ".fasta")
all_files

### read fasta
fasta_file = read.fasta(file = "8_loci_444_spp.fasta", 
                        forceDNACoding = TRUE)

write.dna(fasta_file, 
          file = "your_output_file.phy", 
          format = "interleaved"
          )
