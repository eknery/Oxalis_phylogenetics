### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### chose input directory
dir_input = "3_aligned_sequences/"

### choose output directory
dir_out = "3_aligned_sequences/"

### list locus names
all_loci = list.files(path = paste0(dir_input), pattern = ".fasta")

### looping over loci
for(locus_name in all_loci){
  ### load alignment
  one_locus = read.fasta(paste0( dir_input, locus_name))
  ### species names e number of sites
  spp_names = names(one_locus)
  ### fixed species names
  fixed_spp_names = gsub("_R_", "", spp_names)
  ### replace names with "_R_"
  names(one_locus) = fixed_spp_names
  ### export
  write.fasta(
    sequences = one_locus, 
    as.string = F, 
    names = fixed_spp_names,
    file.out = paste0(dir_out, locus_name),
    nbchar = 100
  )
  ### check!
  print(paste0("Name fixed: ", locus_name)) 
}

