### libraries
if(!require("data.table")) install.packages("data.table"); require("data.table")
if(!require("ape")) install.packages("ape"); library("ape")

### all file names
file_names = list.files("0_accessions")

### list of all files
file_list = list()

### read all files
for(i in 1:length(file_names)){
  file_list[[i]] = read.csv(
    paste0("0_accessions","/",file_names[i]),
    sep = ",",
    h= T
  )
}

### all taxon names
all_tx_names = c()
for(i in 1:length(file_names)){
  all_tx_names = unique(c(all_tx_names, file_list[[i]]$taxon))
}
all_tx_names = sort(all_tx_names)


### joining into a single dt
all_acc = rbindlist(file_list, fill = T)

acc_nums = all_acc$ITS[!is.na(all_acc$ITS)]

x = read.GenBank( 
  access.nb = acc_nums,
  species.names = T
  )

write.table(
  acc_nums, 
  "acc_nums.txt",
  row.names = F,
  quote= F
  )
