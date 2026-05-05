### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### chose input directory
dir_input = "4_trimmed_sequences/"
### choose output directory
dir_out = "5_concatenated_sequences/"
### list locus names
all_loci = list.files(path = paste0(dir_input), pattern = ".fasta")

### select?
locus_filter = "none"
if(locus_filter == "plastid"){
  all_loci = all_loci[!all_loci %in% c("ITS.fasta","ncpGS.fasta")]
}
if(locus_filter == "nuclear"){
  all_loci = all_loci[all_loci %in% c("ITS.fasta","ncpGS.fasta")]
}

### vector to store results
all_partitions = c()
### initial position
ni = 1
for(i in 1:length(all_loci) ){
  ### locus name
  locus_name = all_loci[i]
  ### load alignment
  one_locus = read.fasta(paste0( dir_input, locus_name))
  ### locus length
  locus_length = length(one_locus[[1]])
  ### final position
  nf = (ni - 1) + locus_length
  ### locus partition
  locus_partition = paste0(ni," - ", nf)
  ### update partition matrix
  all_partitions = rbind(all_partitions, c(locus_name, locus_length, locus_partition))
  ### update initial position
  ni = ni + locus_length
}
### into dataframe
all_partitions = as.data.frame(all_partitions)
### column names
colnames(all_partitions) = c("locus", "length", "partition")
### simplify names
all_partitions$locus = gsub(".fasta", "", all_partitions$locus)

### write table
write.table(
  all_partitions,
  paste0(dir_out, locus_filter,"_partition.csv"),
  sep = ",",
  quote = F,
  row.names = F
)

