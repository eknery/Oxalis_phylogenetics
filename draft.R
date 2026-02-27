file_names = list.files("0_accessions")

file_list = list()

for(i in 1:length(file_names))
file_list[[i]] = read.csv(
  paste0("0_accessions","/",file_names[i]),
  sep = ",",
  h= T
  )
