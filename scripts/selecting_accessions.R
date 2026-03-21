### libraries
if(!require("data.table")) install.packages("data.table"); require("data.table")
if(!require("dplyr")) install.packages("dplyr"); library("dplyr")

### all file names
file_names = list.files("0_accessions")

### list of all files
file_list = list()

### read all files
for(i in 1:length(file_names)){
  ### reading file
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

# hasError = c()
# for(i in 1:length(file_names)){
#   hasError = c(hasError, sum(file_list[[i]]$taxon %in% "Oxalis_rhombeo-ovata") )
# }
# names(hasError) = file_names
# hasError

### all column names
all_col_names = c()
for(i in 1:length(file_names)){
  all_col_names = unique( c(all_col_names, colnames(file_list[[i]]) ))
}
all_col_names 

### all marker names
all_marker_names = all_col_names[!all_col_names %in% c("taxon", "specimen")]

### joining into a single dt
all_acc = rbindlist(file_list, fill = T)

### number of entries per marker
n_entry = c()
for(marker_name in all_marker_names){
  n_entry = c(n_entry, sum(!is.na(all_acc[[marker_name]])) )
}
names(n_entry) = all_marker_names

### key marker = maximum number of entries
key_marker_name = names(which.max(n_entry))

### counting sequenced markers and key marker per specimen
marker_index = which(colnames(all_acc) %in% all_marker_names)
all_acc$Nmarker = rowSums(!is.na( all_acc[,..marker_index]))
all_acc$Kmarker = !is.na(all_acc[[key_marker_name]])

### selecting best specimen per taxon
best_acc = c()
for(tx_name in all_tx_names){
  ### select all specimens under a taxon
  tx_acc = all_acc[all_acc$taxon %in% tx_name,]
  ### mark duplicated specimens
  tx_acc$duplicate = duplicated(tx_acc$specimen) | duplicated(tx_acc$specimen, fromLast = TRUE)
  ### if duplicated specimens, join them!
  if( sum(tx_acc$duplicate) != 0 ){
    ### only entries with duplicated specimen
    dup_acc = tx_acc[tx_acc$duplicate == TRUE,]
    ### reference entry
    ref_acc = dup_acc[which.max(dup_acc$Nmarker),]
    ### NA index in reference entry
    na_index = which(is.na(ref_acc))
    ### other entries
    oth_acc = dup_acc[!which.max(dup_acc$Nmarker),]
    ### replace NA for values from other entries
    for(i in na_index){
      for(n in 1:nrow(oth_acc)){
        if(!is.na(oth_acc[n,][[i]]) )
          ref_acc[[i]] = oth_acc[n,][[i]]
      }
    }
    ### only single specimens
    sin_acc = tx_acc[tx_acc$duplicate == FALSE,]
    ### combining single and joined specimens
    tx_acc = rbind(sin_acc, ref_acc)
  }
  ### select entry with maximum number of markers
  max_acc = tx_acc[which.max(tx_acc$Nmarker),]
  ### if more than one entry with maximum number
  if(nrow(max_acc) > 1){
    max_acc = max_acc$Kmarker == TRUE
  }
  ### add selected specimen
  best_acc = rbind(best_acc, max_acc)
}

### export best accession for each taxon
write.table(
  best_acc, 
  paste0("1_best_accessions/","best_acc.csv"),
  sep = ",",
  row.names = F,
  quote= F
  )
