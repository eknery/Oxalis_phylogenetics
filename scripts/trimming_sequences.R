### library
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### chose input directory
dir_input = "3_aligned_sequences/"

### choose output directory
dir_out = "4_trimmed_sequences/"

### list locus names
all_loci = list.files(path = paste0(dir_input), pattern = ".fasta")

### % of missing data to remove a site
threshold = 0.9 #original 0.5

### trimming loci in loop
for(i in 1:length(all_loci) ){
  ### name of one locus
  locus_name = all_loci[i]
  tryCatch({
  ### load alignment
  one_locus = read.fasta(paste0( dir_input, locus_name))
  ### species names e number of sites
  spp_names = names(one_locus)
  n_sites = length(one_locus[[1]])
  ### converting to a matrix
  mtx_locus = matrix(unlist(one_locus), ncol = n_sites, byrow = T)
  ### trimming
  trim = del.colgapsonly(
    x = mtx_locus, 
    threshold = threshold,
    freq.only = FALSE
    )
  ### convert back to matrix
  mtx_trim = as.matrix(as.character(trim))
  ### check if there sites left
  if(ncol(mtx_trim) > 0){
    ### convert back to list
    list_trim = list()
    for(i in 1:nrow(mtx_trim)){
      list_trim[[i]] = paste0(mtx_trim[i,], collapse = "")
    }
    names(list_trim) = spp_names
    ### export
    write.fasta(
      sequences = list_trim, 
      as.string = T, 
      names = spp_names,
      file.out = paste0(dir_out, locus_name),
      nbchar = 100
    )
    ### check!
    print(paste0("Trimming done: ", locus_name)) 
  }
  },
  error = function(e) {
    print(paste0("Skipping: ", locus_name))
    return(NULL)  # Return NULL to indicate failure
  })
}
