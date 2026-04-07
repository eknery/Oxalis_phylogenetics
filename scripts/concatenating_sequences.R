### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### choose directory with sequences
dir_input = "4_trimmed_sequences/" 

### choose directory with sequences
dir_out = "5_concatenated_sequences/"

### list file names
all_loci = list.files(path = paste0(dir_input), pattern = ".fasta")

### plastid only?
plastid = T
if(plastid){
  all_loci = all_loci[!all_loci %in% "ITS.fasta"]
}
### getting all species names across loci
all_spp_names = c()

### loop over loci
for(i in 1:length(all_loci)){
  locus_name = all_loci[i]
  one_locus = read.fasta(paste0(dir_input, locus_name))  
  spp_names = names(one_locus)
  all_spp_names = sort(unique(c(all_spp_names, spp_names)))
}

############################### COMPLETING SEQUENCES ###########################

### list with sequences completed for missing species
complete_loci = list()

### complete sequences with missing species 
for(i in 1:length(all_loci) ){
    ### one locus name
    locus_name = all_loci[i]
    ### pick one locus
    one_locus = read.fasta(paste0(dir_input, locus_name))  
    ### number of sites
    n_sites = length(one_locus[[1]])
    ### species sampled
    sp_names = names(one_locus)
    ### check missing species
    missing_spp = all_spp_names[!all_spp_names %in% sp_names]
    ### including missing species
    if(length(missing_spp) > 0) {
      for(one_missing in missing_spp){
        one_locus[[one_missing]] = rep("-", n_sites)
      }
    }
    ### ordering species
    one_locus = one_locus[all_spp_names]
    ### add locus to final list
    complete_loci[[i]] = one_locus
}

############################### CONCATENATING SEQUENCES ########################

### initial locus
conc_loci = complete_loci[[1]]

### add following loci
for(i in 2:length(complete_loci)){
  adding_locus = complete_loci[[i]]
  for(sp_name in all_spp_names){
    conc_loci[[sp_name]] = c(conc_loci[[sp_name]], adding_locus[[sp_name]])
  }
}

### number of loci
n_loci = length(all_loci)

### number of species
n_spp = length(all_spp_names)

### export
write.fasta(
  sequences = conc_loci, 
  as.string = F, 
  names = all_spp_names,
  file.out = paste0(dir_out, n_loci,"_loci_", n_spp,"_spp.fasta"),
  nbchar = 1000
)
